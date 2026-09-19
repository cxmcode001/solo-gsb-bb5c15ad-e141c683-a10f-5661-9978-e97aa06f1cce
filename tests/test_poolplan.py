"""分池规划测试：pooling 纯算法 + POST /api/batches/{id}/pool-plan 端到端。"""

from __future__ import annotations

import os

import pytest


# ---------- 纯算法 ----------

from app.pooling import CapacityInsufficient, plan_pools


def test_chain_is_conflict_free_when_colorable():
    # 0-1-2 链在 2 池容量 2 下可 0 冲突
    plan = plan_pools(
        ["A", "B", "C"],
        [(0, 1, "high"), (1, 2, "high")],
        pool_count=2,
        pool_capacity=2,
    )
    assert plan.total_conflicts == 0
    assert plan.conflict_edge_count == 2
    assigned = sorted(m for a in plan.assignments for m in a.members)
    assert assigned == [0, 1, 2]
    assert all(a.conflicts == [] for a in plan.assignments)


def test_triangle_forces_one_conflict_and_is_stable():
    edges = [(0, 1, "high"), (1, 2, "medium"), (0, 2, "high")]
    plan = plan_pools(["A", "B", "C"], edges, pool_count=2, pool_capacity=2)
    # 三角形在 2 池中至少 1 条同池边
    assert plan.total_conflicts == 1
    only = plan.assignments[0].conflicts + plan.assignments[1].conflicts
    assert len(only) == 1
    assert (only[0].position_a, only[0].position_b) == (0, 2)
    assert only[0].risk_level == "high"

    # 输入边顺序不同不影响结果
    plan_rev = plan_pools(
        ["A", "B", "C"], list(reversed(edges)), pool_count=2, pool_capacity=2
    )
    assert [a.members for a in plan_rev.assignments] == [
        a.members for a in plan.assignments
    ]
    assert plan_rev.total_conflicts == plan.total_conflicts


def test_vertex_pick_prefers_saturation_then_degree_then_position():
    # 无冲突时：首选取位置 0，并按“新增冲突、池内人数、池号”在两池间轮询
    plan = plan_pools(
        ["A", "B", "C", "D"], [], pool_count=2, pool_capacity=2
    )
    assert plan.assignments[0].members == [0, 2]
    assert plan.assignments[1].members == [1, 3]
    assert plan.total_conflicts == 0


def test_single_pool_with_capacity_collects_all_edges():
    edges = [(0, 1, "high"), (1, 2, "high")]
    plan = plan_pools(["A", "B", "C"], edges, pool_count=1, pool_capacity=3)
    assert len(plan.assignments) == 1
    assert plan.assignments[0].members == [0, 1, 2]
    assert plan.total_conflicts == 2
    assert plan.conflict_edge_count == 2


def test_empty_pools_are_returned():
    plan = plan_pools(["A"], [], pool_count=3, pool_capacity=2)
    assert [a.pool_index for a in plan.assignments] == [0, 1, 2]
    assert plan.assignments[0].members == [0]
    assert plan.assignments[1].members == []
    assert plan.assignments[2].members == []


def test_capacity_insufficient_raises():
    with pytest.raises(CapacityInsufficient) as exc:
        plan_pools(["A", "B", "C"], [], pool_count=2, pool_capacity=1)
    assert exc.value.primer_count == 3
    assert exc.value.total_capacity == 2


def test_capacity_equal_to_primer_count_must_assign_all():
    # 完全图 K3，2 池容量之和恰好 3：仍须全部分配
    plan = plan_pools(
        ["A", "B", "C"],
        [(0, 1, "high"), (0, 2, "high"), (1, 2, "high")],
        pool_count=2,
        pool_capacity=2,
    )
    assigned = sorted(m for a in plan.assignments for m in a.members)
    assert assigned == [0, 1, 2]
    for a in plan.assignments:
        assert len(a.members) <= 2


def test_repeated_calls_are_identical():
    edges = [(0, 3, "high"), (1, 2, "medium"), (2, 3, "high"), (0, 2, "high")]
    first = plan_pools(["P0", "P1", "P2", "P3"], edges, pool_count=3, pool_capacity=2)
    second = plan_pools(["P0", "P1", "P2", "P3"], edges, pool_count=3, pool_capacity=2)
    assert [a.members for a in first.assignments] == [
        a.members for a in second.assignments
    ]
    assert first.total_conflicts == second.total_conflicts


# ---------- API ----------

@pytest.fixture(scope="module")
def client(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("pool_db")
    os.environ["PRIMER_DB_PATH"] = str(tmp_dir / "test.db")

    from app import db, main

    db.init_db()
    from fastapi.testclient import TestClient

    with TestClient(main.app) as c:
        yield c


POOL_BATCH_PAYLOAD = {
    # F=0 R=1 X=2 Q=3 N1=4 N2=5；高风险异二聚体：F-R、F-N1、R-N1（三角形）、X-N2
    "primers": [
        {"name": "F", "sequence": "GGGGGAAAAA"},
        {"name": "R", "sequence": "GGGGGTTTTT"},
        {"name": "X", "sequence": "GGCCCCGG"},
        {"name": "Q", "sequence": "TTGGGGTT"},
        {"name": "N1", "sequence": "AAAATTTT"},
        {"name": "N2", "sequence": "CCCCGGGG"},
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3,
}


@pytest.fixture(scope="module")
def batch_id(client):
    resp = client.post("/api/batches", json=POOL_BATCH_PAYLOAD)
    assert resp.status_code == 201, resp.text
    return resp.json()["id"]


def test_pool_plan_high_only_assigns_all_with_min_conflict(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/pool-plan",
        json={
            "pool_count": 2,
            "pool_capacity": 3,
            "isolate_risk_levels": ["high"],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["batch_id"] == batch_id
    assert data["pool_count"] == 2
    assert data["pool_capacity"] == 3
    assert data["isolate_risk_levels"] == ["high"]
    assert data["primer_count"] == 6
    # 4 条 high 异二聚体边；三角形导致恰有 1 条同池冲突
    assert data["conflict_edge_count"] == 4
    assert data["total_conflicts"] == 1

    all_members = [m for pool in data["pools"] for m in pool["members"]]
    assert sorted(all_members) == ["F", "N1", "N2", "Q", "R", "X"]
    assert all(len(p["members"]) <= 3 for p in data["pools"])

    # 确定性：pool0=[F,N1,N2]（含唯一冲突 F-N1），pool1=[R,X,Q]
    assert data["pools"][0]["members"] == ["F", "N1", "N2"]
    assert data["pools"][0]["member_positions"] == [0, 4, 5]
    pool0_edges = data["pools"][0]["conflicts"]
    assert len(pool0_edges) == 1
    assert pool0_edges[0]["name_a"] == "F"
    assert pool0_edges[0]["name_b"] == "N1"
    assert pool0_edges[0]["risk_level"] == "high"
    assert data["pools"][1]["members"] == ["R", "X", "Q"]
    assert data["pools"][1]["conflicts"] == []

    # 各池 conflicts 数之和等于 total_conflicts
    assert sum(len(p["conflicts"]) for p in data["pools"]) == 1


def test_pool_plan_is_stable_across_calls(client, batch_id):
    body = {
        "pool_count": 3,
        "pool_capacity": 2,
        "isolate_risk_levels": ["high"],
    }
    first = client.post(f"/api/batches/{batch_id}/pool-plan", json=body).json()
    second = client.post(f"/api/batches/{batch_id}/pool-plan", json=body).json()
    assert first["pools"] == second["pools"]
    assert first["total_conflicts"] == second["total_conflicts"]


def test_pool_plan_404_for_unknown_batch(client):
    resp = client.post(
        "/api/batches/no-such-batch/pool-plan",
        json={
            "pool_count": 2,
            "pool_capacity": 3,
            "isolate_risk_levels": ["high"],
        },
    )
    assert resp.status_code == 404
    assert "不存在" in resp.json()["detail"]


def test_pool_plan_422_when_total_capacity_insufficient(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/pool-plan",
        json={
            "pool_count": 2,
            "pool_capacity": 2,  # 总容量 4 < 6
            "isolate_risk_levels": ["high"],
        },
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert "4" in detail and "6" in detail


def test_pool_plan_422_for_bad_params(client, batch_id):
    base = "/api/batches/" + batch_id + "/pool-plan"
    # pool_count 必须 >=1
    assert client.post(
        base,
        json={"pool_count": 0, "pool_capacity": 3, "isolate_risk_levels": ["high"]},
    ).status_code == 422
    # pool_capacity 必须 >=1
    assert client.post(
        base,
        json={"pool_count": 2, "pool_capacity": 0, "isolate_risk_levels": ["high"]},
    ).status_code == 422
    # 空隔离等级列表
    assert client.post(
        base,
        json={"pool_count": 2, "pool_capacity": 3, "isolate_risk_levels": []},
    ).status_code == 422
    # 非法风险等级
    assert client.post(
        base,
        json={"pool_count": 2, "pool_capacity": 3, "isolate_risk_levels": ["critical"]},
    ).status_code == 422
    # 缺字段
    assert client.post(
        base,
        json={"pool_count": 2, "pool_capacity": 3},
    ).status_code == 422


def test_pool_plan_dedupes_risk_levels(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/pool-plan",
        json={
            "pool_count": 2,
            "pool_capacity": 3,
            "isolate_risk_levels": ["HIGH", "high", "High"],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["isolate_risk_levels"] == ["high"]
    assert data["conflict_edge_count"] == 4


def test_pool_plan_medium_level_adds_edges(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/pool-plan",
        json={
            "pool_count": 3,
            "pool_capacity": 2,
            "isolate_risk_levels": ["high", "medium"],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    # high 4 条 + medium 6 条 = 10 条冲突边
    assert data["conflict_edge_count"] == 10
    all_members = [m for pool in data["pools"] for m in pool["members"]]
    assert len(all_members) == 6
    assert len(set(all_members)) == 6
    for edge in (e for p in data["pools"] for e in p["conflicts"]):
        assert edge["risk_level"] in {"high", "medium"}


def test_pool_plan_low_only_and_empty_pools(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/pool-plan",
        json={
            "pool_count": 5,
            "pool_capacity": 2,
            "isolate_risk_levels": ["low"],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    # low 异二聚体 5 条；5 个池全部返回（含空池）
    assert len(data["pools"]) == 5
    assert data["conflict_edge_count"] == 5
    assert sum(len(p["members"]) for p in data["pools"]) == 6
    assert [p["pool_index"] for p in data["pools"]] == [0, 1, 2, 3, 4]


def test_pool_plan_ignores_self_dimers(client):
    # AAAATTTT 的自二聚体高风险，但批次内没有异二聚体 -> 冲突图无边
    created = client.post(
        "/api/batches",
        json={
            "primers": [{"name": "solo", "sequence": "AAAATTTT"}],
            "continuous_threshold": 3,
            "three_prime_threshold": 3,
        },
    ).json()
    resp = client.post(
        f"/api/batches/{created['id']}/pool-plan",
        json={
            "pool_count": 2,
            "pool_capacity": 1,
            "isolate_risk_levels": ["high", "medium", "low"],
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["conflict_edge_count"] == 0
    assert data["total_conflicts"] == 0
    assert data["pools"][0]["members"] == ["solo"]


def test_pool_plan_two_primers_separated_then_forced(client):
    # F/R 为典型 high 3′端异二聚体
    bid = client.post(
        "/api/batches",
        json={
            "primers": [
                {"name": "F", "sequence": "GGGGGAAAAA"},
                {"name": "R", "sequence": "GGGGGTTTTT"},
            ],
            "continuous_threshold": 4,
            "three_prime_threshold": 3,
        },
    ).json()["id"]

    sep = client.post(
        f"/api/batches/{bid}/pool-plan",
        json={
            "pool_count": 2,
            "pool_capacity": 1,
            "isolate_risk_levels": ["high"],
        },
    ).json()
    assert sep["conflict_edge_count"] == 1
    assert sep["total_conflicts"] == 0
    assert sep["pools"][0]["members"] == ["F"]
    assert sep["pools"][1]["members"] == ["R"]

    forced = client.post(
        f"/api/batches/{bid}/pool-plan",
        json={
            "pool_count": 1,
            "pool_capacity": 2,
            "isolate_risk_levels": ["high"],
        },
    ).json()
    assert forced["total_conflicts"] == 1
    assert len(forced["pools"][0]["conflicts"]) == 1
