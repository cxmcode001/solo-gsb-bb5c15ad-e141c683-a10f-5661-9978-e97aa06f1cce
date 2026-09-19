"""冲突图分池规划。

输入异二聚体筛查结果（忽略自二聚体），按需隔离的风险等级把冲突引物对
连成无向冲突图，再用确定性贪心（类 DSATUR）把引物分配到固定数量、固定
容量的池子里：

* 下一个引物的选择键（取最大）：已见邻池数（饱和度）、冲突度（图中的
  总度数）、批次位置（提交顺序，靠前者优先）；
* 放置池子的选择键（取最小）：新增冲突数、池内人数、池号。

选择与放置的全部同分都由批次位置/池号兜底为全序，故相同输入结果稳定、
可复现。仅当 ``pool_count * pool_capacity < 引物数`` 时才会失败，此时
抛出 ``CapacityInsufficient``；容量足够时全部引物必被分配。
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ConflictEdge:
    """冲突图中的一条无向边（一对需隔离的异二聚体）。"""

    position_a: int
    position_b: int
    name_a: str
    name_b: str
    risk_level: str


@dataclass(frozen=True)
class PoolConflictEdge:
    """最终被放进同一个池子的冲突边。"""

    position_a: int
    position_b: int
    name_a: str
    name_b: str
    risk_level: str


@dataclass(frozen=True)
class PoolAssignment:
    pool_index: int
    members: list[int]
    conflicts: list[PoolConflictEdge] = field(default_factory=list)


@dataclass(frozen=True)
class PoolPlan:
    assignments: list[PoolAssignment]
    total_conflicts: int
    conflict_edge_count: int


class CapacityInsufficient(Exception):
    """池子总容量小于引物数量，无法完成全部分配。"""

    def __init__(self, primer_count: int, total_capacity: int) -> None:
        super().__init__(
            f"池子总容量 {total_capacity} 不足以容纳 {primer_count} 个引物"
        )
        self.primer_count = primer_count
        self.total_capacity = total_capacity


def build_conflict_graph(
    positions: list[int],
    edges: list[tuple[int, int, str]],
) -> dict[int, set[int]]:
    """以批次位置为顶点构建无向冲突图。

    ``edges`` 为 ``(position_a, position_b, risk_level)``，调用方需保证
    两个位置都存在于 ``positions`` 且为异二聚体对。
    """
    graph: dict[int, set[int]] = {pos: set() for pos in positions}
    for a, b, _risk in edges:
        graph[a].add(b)
        graph[b].add(a)
    return graph


def plan_pools(
    names: list[str],
    hetero_conflicts: list[tuple[int, int, str]],
    *,
    pool_count: int,
    pool_capacity: int,
) -> PoolPlan:
    """按确定性贪心把 ``names``（索引即批次位置）分配到各池。

    ``hetero_conflicts`` 为 ``(位置 a, 位置 b, 风险等级)``，仅包含命中
    需隔离风险等级的异二聚体对（自二聚体已忽略）。
    """
    n = len(names)
    if pool_count * pool_capacity < n:
        raise CapacityInsufficient(n, pool_count * pool_capacity)

    positions = list(range(n))
    graph = build_conflict_graph(positions, hetero_conflicts)

    # 边元信息（同一条边只存一次，位置小者在前）
    edge_info: dict[tuple[int, int], str] = {}
    for a, b, risk in hetero_conflicts:
        key = (a, b) if a < b else (b, a)
        edge_info.setdefault(key, risk)

    # 已放置顶点 -> 池号
    placed: dict[int, int] = {}
    # 池号 -> 成员（位置）
    pools: list[list[int]] = [[] for _ in range(pool_count)]
    # 顶点 -> 邻接顶点已落入的池号集合（饱和度依据）
    neighbor_pools: dict[int, set[int]] = {pos: set() for pos in positions}

    for _ in range(n):
        # 选择下一个引物：已见邻池数降序、冲突度（总度数）降序、批次位置升序
        vertex = max(
            (pos for pos in positions if pos not in placed),
            key=lambda pos: (
                len(neighbor_pools[pos]),
                len(graph[pos]),
                -pos,
            ),
        )

        # 选择池子：新增冲突数升序、池内人数升序、池号升序；容量满的池跳过
        def pool_key(pool: int) -> tuple[int, int, int]:
            added = sum(
                1 for nb in graph[vertex] if placed.get(nb) == pool
            )
            return (added, len(pools[pool]), pool)

        target = min(
            (p for p in range(pool_count) if len(pools[p]) < pool_capacity),
            key=pool_key,
        )

        placed[vertex] = target
        pools[target].append(vertex)
        for nb in graph[vertex]:
            if nb not in placed:
                neighbor_pools[nb].add(target)

    assignments: list[PoolAssignment] = []
    total_conflicts = 0
    for pool, members in enumerate(pools):
        member_set = set(members)
        pool_edges: list[PoolConflictEdge] = []
        for (a, b), risk in edge_info.items():
            if a in member_set and b in member_set:
                pool_edges.append(
                    PoolConflictEdge(
                        position_a=a,
                        position_b=b,
                        name_a=names[a],
                        name_b=names[b],
                        risk_level=risk,
                    )
                )
        pool_edges.sort(key=lambda e: (e.position_a, e.position_b))
        members_sorted = sorted(members)
        total_conflicts += len(pool_edges)
        assignments.append(
            PoolAssignment(
                pool_index=pool,
                members=members_sorted,
                conflicts=pool_edges,
            )
        )

    return PoolPlan(
        assignments=assignments,
        total_conflicts=total_conflicts,
        conflict_edge_count=len(edge_info),
    )
