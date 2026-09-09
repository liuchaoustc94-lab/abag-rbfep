"""Mutation-route classification tests (user policy table 2026-08-26)."""

from dataclasses import dataclass

from abag_rbfe.routing import (
    ROUTE_BASELINE,
    ROUTE_DSSB,
    ROUTE_ENDPOINT_ENSEMBLE,
    ROUTE_PROTONATION,
    ROUTE_RID_HYDRATION,
    classify_mutation_route,
)


@dataclass(frozen=True)
class _Site:
    wt: str
    mut: str
    resseq: int = 50


def _route(sites, charge_conserving=True, loop_residues=None):
    return classify_mutation_route(
        tuple(sites),
        charge_conserving=charge_conserving,
        loop_residues=loop_residues,
    )


def test_ordinary_ala_scan_is_baseline() -> None:
    route = _route([_Site("S", "A")])
    assert route.primary == ROUTE_BASELINE
    assert route.tags == ()


def test_aromatic_deletion_routes_to_rid_hydration() -> None:
    for wt in ("Y", "W", "F", "M"):
        route = _route([_Site(wt, "A")])
        assert route.primary == ROUTE_RID_HYDRATION, wt
        assert "basin_diversity" in route.tags


def test_partial_size_reduction_is_not_cavity_deletion() -> None:
    # Y -> F keeps the aromatic ring: not a cavity-forming deletion.
    route = _route([_Site("Y", "F")])
    assert route.primary == ROUTE_BASELINE


def test_pro_gly_and_insertion_route_to_endpoint_ensemble() -> None:
    assert _route([_Site("P", "V")]).primary == ROUTE_ENDPOINT_ENSEMBLE
    assert _route([_Site("G", "W")]).primary == ROUTE_ENDPOINT_ENSEMBLE
    # plain insertion (small -> large, no Pro/Gly)
    assert _route([_Site("A", "M")]).primary == ROUTE_ENDPOINT_ENSEMBLE


def test_loop_annotation_forces_endpoint_ensemble() -> None:
    route = _route([_Site("S", "T", 30)], loop_residues=frozenset({30}))
    assert route.primary == ROUTE_ENDPOINT_ENSEMBLE
    assert "flexible_loop" in route.tags


def test_charge_change_routes_to_dssb_with_priority() -> None:
    # charge change + aromatic deletion: DSSB wins as primary
    route = _route([_Site("Y", "D")], charge_conserving=False)
    assert route.primary == ROUTE_DSSB


def test_titratable_neighborhood_routes_to_protonation() -> None:
    route = _route([_Site("H", "Q")])
    assert route.primary == ROUTE_PROTONATION
    assert "titratable_neighbor" in route.tags


def test_titratable_is_tag_when_higher_priority_applies() -> None:
    # K -> D is charge-changing: DSSB primary, titratable stays a tag
    route = _route([_Site("K", "D")], charge_conserving=False)
    assert route.primary == ROUTE_DSSB
    assert "titratable_neighbor" in route.tags


def test_deletion_of_small_residue_stays_baseline() -> None:
    route = _route([_Site("V", "A")])
    assert route.primary == ROUTE_BASELINE
