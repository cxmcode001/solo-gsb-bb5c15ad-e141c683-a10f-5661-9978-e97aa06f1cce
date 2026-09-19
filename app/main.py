"""FastAPI 入口：批次提交、批次结果查询、批次列表。"""

from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import RedirectResponse

from . import config, db
from .schemas import (
    BatchCreate,
    BatchList,
    BatchResult,
    BatchSummary,
    PairOut,
    PoolConflictEdgeOut,
    PoolOut,
    PoolPlanRequest,
    PoolPlanResult,
    PrimerOut,
)
from .services import build_pool_plan, run_batch


@asynccontextmanager
async def lifespan(_app: FastAPI):
    # 依赖与数据库初始化纳入容器启动流程
    db.init_db(reset=config.RESET_ON_START)
    yield


app = FastAPI(
    title="PCR 引物二聚体筛查 API",
    description=(
        "提交一批带名称的 IUPAC DNA 引物，枚举自二聚体与两两异二聚体的全部"
        "错位对齐，按 IUPAC 碱基集合判断互补，返回每对最危险对齐、连续互补长度、"
        "涉及 3′端的连续长度与风险原因。"
    ),
    version="1.0.0",
    lifespan=lifespan,
)


@app.get("/", include_in_schema=False)
def root() -> RedirectResponse:
    return RedirectResponse(url="/docs")


@app.get("/health", tags=["meta"], summary="健康检查")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post(
    "/api/batches",
    status_code=201,
    tags=["batches"],
    summary="提交引物批次并同步完成筛查",
    response_model=BatchResult,
)
def create_batch(payload: BatchCreate) -> BatchResult:
    if len(payload.primers) > config.MAX_PRIMERS_PER_BATCH:
        raise HTTPException(
            status_code=422,
            detail=(
                f"引物数量 {len(payload.primers)} 超过上限 "
                f"{config.MAX_PRIMERS_PER_BATCH}"
            ),
        )

    seen: set[str] = set()
    duplicates: list[str] = []
    primers: list[tuple[str, str]] = []
    for p in payload.primers:
        key = p.name
        if key in seen:
            duplicates.append(key)
        seen.add(key)
        primers.append((p.name, p.sequence))
    if duplicates:
        raise HTTPException(
            status_code=422,
            detail=f"批次内引物名称重复: {', '.join(duplicates)}",
        )

    batch_id = run_batch(
        primers=primers,
        continuous_threshold=payload.continuous_threshold,
        three_prime_threshold=payload.three_prime_threshold,
    )
    return _load_result(batch_id, risky_only=False)  # type: ignore[return-value]


@app.get(
    "/api/batches",
    tags=["batches"],
    summary="列出历史批次（按创建时间倒序）",
    response_model=BatchList,
)
def list_batches(
    limit: int = Query(50, ge=1, le=200),
    offset: int = Query(0, ge=0),
) -> BatchList:
    with db.get_conn() as conn:
        total = db.count_batches(conn)
        rows = db.list_batches(conn, limit, offset)
    return BatchList(
        total=total,
        limit=limit,
        offset=offset,
        batches=[_summary(r) for r in rows],
    )


@app.get(
    "/api/batches/{batch_id}",
    tags=["batches"],
    summary="查询批次完整筛查结果",
    response_model=BatchResult,
    responses={404: {"description": "批次不存在"}},
)
def get_batch(
    batch_id: str,
    risky_only: bool = Query(
        False, description="为 true 时只返回 high/medium 风险配对"
    ),
) -> BatchResult:
    return _load_result(batch_id, risky_only=risky_only)


@app.post(
    "/api/batches/{batch_id}/pool-plan",
    tags=["batches"],
    summary="基于异二聚体冲突图为已保存批次规划引物分池",
    response_model=PoolPlanResult,
    responses={
        404: {"description": "批次不存在"},
        422: {"description": "池子总容量不足或请求参数非法"},
    },
)
def create_pool_plan(batch_id: str, payload: PoolPlanRequest) -> PoolPlanResult:
    # 先确认批次存在，再给出容量结论，错误语义与文档一致
    with db.get_conn() as conn:
        batch = db.get_batch(conn, batch_id)
        if batch is None:
            raise HTTPException(status_code=404, detail="批次不存在")
        primer_count = batch["primer_count"]

    total_capacity = payload.pool_count * payload.pool_capacity
    if total_capacity < primer_count:
        raise HTTPException(
            status_code=422,
            detail=(
                f"池子总容量 {total_capacity}（{payload.pool_count} 池 × "
                f"{payload.pool_capacity}）不足以容纳 {primer_count} 个引物"
            ),
        )

    plan = build_pool_plan(
        batch_id,
        pool_count=payload.pool_count,
        pool_capacity=payload.pool_capacity,
        isolate_risk_levels=payload.isolate_risk_levels,
    )

    with db.get_conn() as conn:
        primers = db.get_primers(conn, batch_id)
    pos_to_name = {row["position"]: row["name"] for row in primers}

    pools = [
        PoolOut(
            pool_index=assignment.pool_index,
            members=[pos_to_name[pos] for pos in assignment.members],
            member_positions=assignment.members,
            conflicts=[
                PoolConflictEdgeOut(
                    position_a=edge.position_a,
                    position_b=edge.position_b,
                    name_a=edge.name_a,
                    name_b=edge.name_b,
                    risk_level=edge.risk_level,
                )
                for edge in assignment.conflicts
            ],
        )
        for assignment in plan.assignments
    ]

    return PoolPlanResult(
        batch_id=batch_id,
        pool_count=payload.pool_count,
        pool_capacity=payload.pool_capacity,
        isolate_risk_levels=payload.isolate_risk_levels,
        primer_count=primer_count,
        conflict_edge_count=plan.conflict_edge_count,
        total_conflicts=plan.total_conflicts,
        pools=pools,
    )


def _summary(row) -> BatchSummary:
    return BatchSummary(**dict(row))


def _load_result(batch_id: str, *, risky_only: bool) -> BatchResult:
    with db.get_conn() as conn:
        batch = db.get_batch(conn, batch_id)
        if batch is None:
            raise HTTPException(status_code=404, detail="批次不存在")
        primers = db.get_primers(conn, batch_id)
        pairs = db.get_pairs(conn, batch_id, risky_only=risky_only)

    return BatchResult(
        **dict(batch),
        primers=[PrimerOut(**dict(p)) for p in primers],
        pairs=[PairOut(**db.serialize_pair(r)) for r in pairs],
    )
