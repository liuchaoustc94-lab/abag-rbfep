from pathlib import Path

import pytest

from abag_pmx.mutations import (
    build_job_id,
    build_mutation_group,
    load_mutation_groups_from_csv,
    mutation_script_lines,
    parse_mutation_group_tokens,
)


def test_single_and_same_side_double_groups_are_supported(tmp_path: Path) -> None:
    csv_path = tmp_path / "mutations.csv"
    csv_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
                "double_h_pair,H,52,,S,T,antibody",
                "double_h_pair,H,54,,N,Q,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    groups = load_mutation_groups_from_csv(csv_path)
    single = build_mutation_group(groups[0]["mutation_group_id"], groups[0]["sites"], True, False)
    double = build_mutation_group(groups[1]["mutation_group_id"], groups[1]["sites"], True, False)

    assert single.mutation_count == 1
    assert single.min_version == "v1"
    assert double.mutation_count == 2
    assert double.min_version == "v2"
    assert double.entity_side == "antibody"


def test_cross_side_double_is_rejected() -> None:
    sites = parse_mutation_group_tokens("H:Y32F@antibody;A:T52S@antigen")
    with pytest.raises(ValueError, match="same-side"):
        build_mutation_group("cross_side", sites, True, False)


def test_pmx_script_lines_reject_insertion_codes() -> None:
    sites = parse_mutation_group_tokens("H:Y32AF@antibody")
    with pytest.raises(ValueError, match="insertion codes"):
        mutation_script_lines(tuple(sites))


def test_job_id_is_stable_and_readable() -> None:
    group = build_mutation_group("single_h_y32f", parse_mutation_group_tokens("H:Y32F@antibody"), True, False)
    assert build_job_id("Demo_ABAG", group).startswith("demo-abag-antibody-")

