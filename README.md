# PCR 引物二聚体筛查 API

基于 **Python + FastAPI + SQLite** 的引物二聚体（primer dimer）批量筛查服务。
调用方提交一批带名称的 IUPAC DNA 引物及两个阈值（连续互补阈值、3′端连续互补阈值），
服务枚举每对引物**自二聚体**与**两两异二聚体的全部错位对齐**，按 IUPAC 简并碱基
**集合语义**判断碱基互补，返回每对的最危险对齐、最长连续互补长度、涉及任一 3′端的
连续长度及中文风险原因；同分时按固定偏移规则取唯一结果。批次与结果持久化到
SQLite，支持事后查询。

---

## 一条命令启动

```bash
docker compose up --build -d
```

启动后服务监听宿主机 **8000 端口**（`docker-compose.yml` 中 `8000:8000`
端口映射），在本机浏览器或命令行直接访问 `http://localhost:8000/...` 即可，
无需额外端口或反向代理。

访问入口：

| 用途 | 地址 |
| --- | --- |
| 交互式 API 文档（Swagger UI） | http://localhost:8000/docs |
| ReDoc 文档 | http://localhost:8000/redoc |
| 健康检查 | http://localhost:8000/health |
| 批次接口 | http://localhost:8000/api/batches |

> 监听端口 **8000**。SQLite 数据库文件保存在命名卷 `primer-data` 内
> （容器路径 `/data/primer_dimer.db`）。镜像构建期安装依赖、容器启动期自动建库
> （`docker-entrypoint.sh` → `python -m app.cli initdb` → `uvicorn`）。

停止与查看日志：

```bash
docker compose logs -f
docker compose down            # 保留数据卷
docker compose down -v         # 同时删除数据库卷
```

---

## 必要配置（环境变量，均有默认值，不配也能直接启动）

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `PRIMER_DB_PATH` | `/data/primer_dimer.db` | SQLite 数据库文件路径 |
| `PRIMER_HOST` | `0.0.0.0` | 监听地址 |
| `PRIMER_PORT` | `8000` | 监听端口（容器内；改宿主机端口请编辑 `docker-compose.yml` 的 `ports`） |
| `PRIMER_MAX_PRIMERS` | `500` | 单批次引物数量上限 |
| `PRIMER_MAX_LENGTH` | `200` | 单条引物长度上限 |
| `PRIMER_DB_RESET` | `0` | 设为 `1` 时启动清空并重建表（仅用于重置数据） |

---

## API 概览

### `POST /api/batches` — 提交批次（同步计算后返回完整结果）

```bash
curl -s -X POST http://localhost:8000/api/batches \
  -H 'Content-Type: application/json' \
  -d '{
    "primers": [
      {"name": "F", "sequence": "GGGGGAAAAA"},
      {"name": "R", "sequence": "GGGGGTTTTT"}
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3
  }'
```

请求字段：

| 字段 | 说明 |
| --- | --- |
| `primers[].name` | 引物名称，非空，批次内唯一 |
| `primers[].sequence` | 5′→3′ IUPAC DNA 序列，接受 `A T C G R Y S W K M B D H V N`（大小写、首尾空白自动规范化；含 `-`、空格等非 IUPAC 字符返回 422） |
| `continuous_threshold` | 连续互补长度阈值（≥1），达到即提示双链结合稳定 |
| `three_prime_threshold` | 3′端连续互补长度阈值（≥1），达到即提示聚合酶延伸风险，风险等级 `high` |

返回（节选）：

```json
{
  "id": "f98cdd27-...",
  "primer_count": 2,
  "pair_count": 3,
  "high_risk_count": 1, "medium_risk_count": 0, "low_risk_count": 2,
  "pairs": [
    {
      "pair_index": 2,
      "name_top": "F", "name_bottom": "R", "is_self_dimer": false,
      "alignments_checked": 19,
      "best_shift": 5,
      "overlap_len": 10,
      "continuous_complementary_length": 5,
      "three_prime_run_length": 5,
      "top3_run_length": 5, "bottom3_run_length": 5,
      "involved_three_prime_ends": ["F:3′", "R:3′"],
      "risk_level": "high",
      "risk_reasons": [
        "3′端连续互补达到阈值，存在聚合酶延伸风险：top 3′端连续互补 5 bp，bottom 3′端连续互补 5 bp",
        "最长连续互补 5 bp 达到阈值 4 bp，双链结合稳定"
      ],
      "alignment": {
        "top":    "GGGGGAAAAA     ",
        "match":  "     |||||     ",
        "bottom": "     TTTTTGGGGG"
      }
    }
  ]
}
```

### `GET /api/batches/{id}` — 批次结果查询

* `GET /api/batches/{id}` 返回该批次全部配对；
* `GET /api/batches/{id}?risky_only=true` 只返回 `high`/`medium` 风险配对；
* 不存在的批次返回 `404`。

### `GET /api/batches?limit=50&offset=0` — 批次列表（按创建时间倒序）

### `POST /api/batches/{id}/pool-plan` — 基于异二聚体冲突图的引物分池

对**已保存批次**做分池规划：忽略全部自二聚体，把风险等级命中
`isolate_risk_levels` 的**异二聚体**对连成无向冲突图，再用确定性贪心
（类 DSATUR）把引物分入指定数量、统一容量的池子。

* **选择下一个引物**（取最大）：已见邻池数（邻接引物已落入的不同池数）
  → 冲突度（冲突图总度数）→ 批次位置（提交顺序靠前者优先）；
* **选择放置池子**（取最小）：新增冲突数（与池内成员的冲突边数）
  → 池内人数 → 池号；
* 只要 `pool_count × pool_capacity ≥ 引物数`，全部引物必被分配；
  同分全部由批次位置/池号兜底，**相同输入结果稳定可复现**。

```bash
curl -s -X POST http://localhost:8000/api/batches/$BID/pool-plan \
  -H 'Content-Type: application/json' \
  -d '{
    "pool_count": 2,
    "pool_capacity": 3,
    "isolate_risk_levels": ["high"]
  }'
```

请求字段：

| 字段 | 说明 |
| --- | --- |
| `pool_count` | 池子数量，≥1（≤1000） |
| `pool_capacity` | 每个池子的容量，≥1（≤100000），各池相同 |
| `isolate_risk_levels` | 需隔离的风险等级列表，取值 `high` / `medium` / `low`；异二聚体命中其中任一等级即连冲突边；大小写不敏感、重复值自动去重，至少提供一个 |

返回：各池成员（名称 + 批次位置）、**同池冲突边**（两端引物的名称、位置与风险等级）
以及**总冲突数**（各池冲突边数之和）；空池同样返回。

```json
{
  "batch_id": "115fdb0d-...",
  "pool_count": 2,
  "pool_capacity": 3,
  "isolate_risk_levels": ["high"],
  "primer_count": 6,
  "conflict_edge_count": 4,
  "total_conflicts": 1,
  "pools": [
    {
      "pool_index": 0,
      "members": ["F", "N1", "N2"],
      "member_positions": [0, 4, 5],
      "conflicts": [
        {"position_a": 0, "position_b": 4, "name_a": "F", "name_b": "N1", "risk_level": "high"}
      ]
    },
    {
      "pool_index": 1,
      "members": ["R", "X", "Q"],
      "member_positions": [1, 2, 3],
      "conflicts": []
    }
  ]
}
```

错误：

* 批次不存在：`404`（`{"detail": "批次不存在"}`）；
* 池子总容量不足（`pool_count × pool_capacity < 引物数`）：`422`，
  响应说明总容量与引物数，例如
  `池子总容量 4（2 池 × 2）不足以容纳 6 个引物`；
* 参数非法（池数/容量小于 1、空或非法的风险等级列表等）：`422`。

> 示例中的 `$BID` 为 `POST /api/batches` 返回的批次 `id`。
> 所有接口均通过 **8000** 端口访问，例如
> http://localhost:8000/api/batches 与 http://localhost:8000/docs 。

## 配对结果字段含义

| 字段 | 含义 |
| --- | --- |
| `is_self_dimer` | 是否为自二聚体（同一引物与自身对齐） |
| `alignments_checked` | 实际枚举的错位数量，恒为 `len(top)+len(bottom)-1` |
| `best_shift` | 最危险对齐的带符号错位（见下） |
| `overlap_len` | 该对齐重叠碱基数 |
| `continuous_complementary_length` | 重叠区内最长连续互补长度 |
| `top3_run_length` / `bottom3_run_length` | 触及对应 3′端的连续互补长度（端未参与则为 0） |
| `three_prime_run_length` | 两端的最大值（`max(top3, bottom3)`） |
| `involved_three_prime_ends` | 3′端连续互补长度 > 0 的端名称 |
| `risk_level` | `high`（3′端达阈值）/ `medium`（仅连续互补达阈值）/ `low` |
| `risk_reasons` | 中文风险原因列表 |
| `alignment` | 等宽三行可视化；中间行 `\|`=互补，`.`=不互补，空格=未重叠 |

## 错位 `shift` 定义与同分规则

两条引物按反向平行对齐，top 固定 5′→3′，bottom 按 3′→5′ 打印。
`shift`（记为 d）是 **bottom 3′端碱基所在的全局列坐标**，top 的 5′端 T0 位于列 0：

下列三行与 API 的 `alignment` 输出格式一致（等宽，空格为未重叠，`|` 为互补）：

```
d = 0（齐平）        d > 0（bottom 3′伸入 top）   d < 0（bottom 3′悬空于左侧）
T0 T1 T2             T0 T1 T2                        T0 T1 T2
|  |  |              .  |  |  |                  |  |  .
   B2 B1 B0             B2 B1 B0              B2 B1 B0
（渲染为连续字符行，例如 d=0 时 top="TTT"、bottom="BBB"）
```

* `d = 0`：bottom 3′端碱基与 top 5′端 T0 同列；
* `d > 0`：bottom 3′端向 top 3′方向伸入 top 内部；
* `d < 0`：bottom 3′端悬空于 top 5′端左侧。

枚举范围为 `-(b-1) <= d <= a-1`，共 `a+b-1` 个对齐，全部检查，不做启发式剪枝。

**最危险对齐**按以下键取最大者：

1. `max(top3_run, bottom3_run)`（聚合酶延伸风险，最重要）；
2. 最长连续互补长度（双链稳定性）；
3. 重叠碱基数。

**同分固定偏移规则**（保证结果唯一、可复现）：风险键完全相同时依次取

1. `|d|` 最小者（悬空端最短、最稳定构象）；
2. `d` 为正者优先，其次 0，最后负值；
3. 再按 `d` 升序兜底（全序，无并列残留）。

## IUPAC 简并碱基互补规则

每个简并码映射为碱基集合（如 `R={A,G}`、`Y={C,T}`、`N={A,C,G,T}`）。
同列两碱基“可能互补”当且仅当一个碱基集合的 Watson-Crick 互补集合与另一个
集合**交集非空**。因此 `R` 与 `Y` 互补、`N` 与任意碱基互补、`S={G,C}` 与
`W={A,T}` 不互补。只要存在至少一种简并实例化方式可形成 WC 配对即判为互补。

---

## 本地开发（不使用 Docker）

```bash
pip install -r requirements.txt
export PRIMER_DB_PATH=./data/primer_dimer.db
python -m app.cli initdb
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试：

```bash
pip install pytest
python -m pytest -q
```

## 项目结构

```
app/
  config.py     # 环境变量配置
  iupac.py      # IUPAC 碱基集合与互补矩阵
  dimer.py      # 错位枚举、连续互补/3′端长度、风险判定、同分规则
  pooling.py    # 异二聚体冲突图构建与确定性贪心分池（pool-plan）
  schemas.py    # Pydantic 模型与输入校验
  db.py         # SQLite 建表、写入、查询
  services.py   # 批次计算编排（N 自二聚体 + N(N-1)/2 异二聚体）与分池加载
  main.py       # FastAPI 路由
  cli.py        # initdb 命令
tests/          # pytest 算法与 API 测试
Dockerfile
docker-compose.yml
docker-entrypoint.sh
requirements.txt
```
