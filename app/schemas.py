"""Pydantic 请求/响应模型与输入校验。"""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator

from . import config
from .iupac import validate_iupac


class PrimerIn(BaseModel):
    name: str = Field(..., description="引物名称（非空，批次内唯一）")
    sequence: str = Field(
        ..., description="5'→3' IUPAC DNA 序列（允许 ATCGRYSWKMBDHVN）"
    )

    @field_validator("name")
    @classmethod
    def _name(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("引物名称不能为空")
        if len(v) > 100:
            raise ValueError("引物名称最长 100 字符")
        return v

    @field_validator("sequence")
    @classmethod
    def _sequence(cls, v: str) -> str:
        seq = validate_iupac(v)
        if len(seq) > config.MAX_PRIMER_LENGTH:
            raise ValueError(
                f"引物长度 {len(seq)} 超过上限 {config.MAX_PRIMER_LENGTH}"
            )
        return seq


class BatchCreate(BaseModel):
    primers: list[PrimerIn] = Field(
        ..., min_length=1, description="带名称的 IUPAC DNA 引物批次"
    )
    continuous_threshold: int = Field(
        ...,
        ge=1,
        le=config.MAX_PRIMER_LENGTH,
        description="连续互补长度阈值（达到即提示结合稳定风险）",
    )
    three_prime_threshold: int = Field(
        ...,
        ge=1,
        le=config.MAX_PRIMER_LENGTH,
        description="3′端连续互补长度阈值（达到即提示延伸风险，风险等级 high）",
    )


class AlignmentOut(BaseModel):
    top: str
    match: str
    bottom: str


class PairOut(BaseModel):
    pair_index: int
    name_top: str
    name_bottom: str
    is_self_dimer: bool
    alignments_checked: int
    best_shift: int = Field(
        description=(
            "最危险对齐的带符号错位：bottom 3′端碱基的全局列坐标（top 5′端 "
            "T0 位于列 0）。0=齐平，>0=bottom 3′端伸入 top，<0=悬空于 top 5′左侧"
        )
    )
    overlap_len: int
    continuous_complementary_length: int
    three_prime_run_length: int
    top3_run_length: int
    bottom3_run_length: int
    involved_three_prime_ends: list[str]
    risk_level: str = Field(description="high / medium / low")
    risk_reasons: list[str]
    alignment: AlignmentOut


class PrimerOut(BaseModel):
    position: int
    name: str
    sequence: str


class BatchSummary(BaseModel):
    id: str
    created_at: str
    continuous_threshold: int
    three_prime_threshold: int
    primer_count: int
    pair_count: int
    high_risk_count: int
    medium_risk_count: int
    low_risk_count: int


class BatchResult(BatchSummary):
    primers: list[PrimerOut]
    pairs: list[PairOut]


class BatchList(BaseModel):
    total: int
    limit: int
    offset: int
    batches: list[BatchSummary]


_RISK_LEVELS = {"high", "medium", "low"}


class PoolPlanRequest(BaseModel):
    pool_count: int = Field(
        ..., ge=1, le=1000, description="池子数量（≥1）"
    )
    pool_capacity: int = Field(
        ..., ge=1, le=100000, description="每个池子的容量（≥1，各池相同）"
    )
    isolate_risk_levels: list[str] = Field(
        ...,
        min_length=1,
        description=(
            "需隔离的风险等级（异二聚体命中其中任一等级即在冲突图中连边，"
            "自二聚体忽略）；取值 high / medium / low，重复值自动去重"
        ),
    )

    @field_validator("isolate_risk_levels")
    @classmethod
    def _isolate_risk_levels(cls, v: list[str]) -> list[str]:
        levels = [s.strip().lower() for s in v]
        bad = sorted({s for s in levels if s not in _RISK_LEVELS})
        if bad:
            raise ValueError(
                "风险等级仅支持 high / medium / low，非法值: " + ", ".join(bad)
            )
        # 去重并保持首次出现顺序，保证相同语义输入结果稳定
        return list(dict.fromkeys(levels))


class PoolConflictEdgeOut(BaseModel):
    position_a: int
    position_b: int
    name_a: str
    name_b: str
    risk_level: str


class PoolOut(BaseModel):
    pool_index: int
    members: list[str] = Field(
        description="池内引物名称（按批次位置升序）"
    )
    member_positions: list[int] = Field(
        description="池内引物在批次中的位置（与 members 一一对应）"
    )
    conflicts: list[PoolConflictEdgeOut] = Field(
        description="同池冲突边（两端引物被分入本池的冲突对）"
    )


class PoolPlanResult(BaseModel):
    batch_id: str
    pool_count: int
    pool_capacity: int
    isolate_risk_levels: list[str]
    primer_count: int
    conflict_edge_count: int = Field(
        description="冲突图中的边数（命中隔离等级的异二聚体对数）"
    )
    total_conflicts: int = Field(
        description="最终同池冲突边总数（各池 conflicts 长度之和）"
    )
    pools: list[PoolOut]


class ErrorOut(BaseModel):
    detail: str
