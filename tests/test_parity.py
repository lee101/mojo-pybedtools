from __future__ import annotations

import pybedtools
import pytest

from mojo_pybedtools import BedTool, Interval, create_interval_from_list


A_TEXT = """chr1\t0\t10\ta\t0\t+
chr1\t20\t30\tb\t0\t-
chr2\t0\t5\tc\t0\t.
"""
B_TEXT = """chr1\t5\t25\tx\t0\t+
chr1\t8\t12\ty\t0\t-
chr2\t3\t9\tz\t0\t.
"""


def ours(text: str) -> BedTool:
    return BedTool(text, from_string=True)


def upstream(text: str):
    return pybedtools.BedTool(text, from_string=True)


@pytest.mark.parametrize("kwargs", [
    {}, {"wa": True}, {"wb": True}, {"wa": True, "wb": True}, {"wo": True},
    {"u": True}, {"v": True}, {"c": True}, {"s": True}, {"S": True},
    {"f": 0.5}, {"F": 0.5}, {"f": 0.5, "F": 0.5, "e": True}, {"f": 0.5, "r": True},
])
def test_intersect_matches_pybedtools(kwargs):
    assert str(ours(A_TEXT).intersect(ours(B_TEXT), **kwargs)) == str(upstream(A_TEXT).intersect(upstream(B_TEXT), **kwargs))


def test_subtract_and_drop_any_overlap_match_pybedtools():
    assert str(ours(A_TEXT).subtract(ours(B_TEXT))) == str(upstream(A_TEXT).subtract(upstream(B_TEXT)))
    assert str(ours(A_TEXT).subtract(ours(B_TEXT), A=True)) == str(upstream(A_TEXT).subtract(upstream(B_TEXT), A=True))


@pytest.mark.parametrize("kwargs", [
    {"f": 0.5}, {"F": 0.5}, {"f": 0.5, "r": True},
    {"f": 0.5, "F": 0.5, "e": True}, {"s": True}, {"S": True},
])
def test_subtract_overlap_filters_match_pybedtools(kwargs):
    assert str(ours(A_TEXT).subtract(ours(B_TEXT), **kwargs)) == str(upstream(A_TEXT).subtract(upstream(B_TEXT), **kwargs))


def test_coverage_deduplicates_overlapping_b_features():
    a = "chr1\t0\t20\ta\t0\t+\n"
    b = "chr1\t2\t12\tx\t0\t+\nchr1\t8\t18\ty\t0\t+\n"
    assert str(ours(a).coverage(ours(b))) == str(upstream(a).coverage(upstream(b)))


def test_coverage_preserves_a_input_order_and_handles_contained_b_features():
    a = "chr1\t200\t210\tb\nchr1\t0\t100\ta\n"
    b = "chr1\t0\t100\tx\nchr1\t50\t60\ty\n"
    assert str(ours(a).coverage(ours(b))) == str(upstream(a).coverage(upstream(b)))


def test_intersect_count_handles_empty_and_nonempty_records():
    a = "chr1\t0\t10\na\t0\t0\nchr1\t20\t30\n"
    b = "chr1\t5\t8\nchr1\t21\t22\n"
    assert str(ours(a).intersect(ours(b), c=True)) == str(upstream(a).intersect(upstream(b), c=True))


def test_merge_matches_pybedtools_for_unsorted_and_bookended_input():
    text = "chr2\t8\t10\nchr1\t8\t12\nchr1\t0\t5\nchr1\t5\t8\nchr1\t20\t22\n"
    assert str(ours(text).merge()) == str(upstream(text).sort().merge())
    assert str(ours(text).merge(d=8)) == str(upstream(text).sort().merge(d=8))


def test_empty_merge_does_not_cross_the_ffi_boundary():
    assert str(ours("").merge()) == ""


def test_materialized_results_remain_valid_inputs_to_kernels():
    a, b = ours(A_TEXT), ours(B_TEXT)
    counted = a.intersect(b, c=True)
    covered = a.coverage(b)
    merged = a.merge()
    assert str(counted.intersect(b, c=True)) == str(upstream(str(counted)).intersect(upstream(B_TEXT), c=True))
    assert str(covered.coverage(b)) == str(upstream(str(covered)).coverage(upstream(B_TEXT)))
    assert str(merged.merge()) == str(upstream(str(merged)).sort().merge())


def test_interval_surface_and_file_roundtrip(tmp_path):
    interval = create_interval_from_list(["chr4", 7, 11, "n", 3, "+"])
    assert isinstance(interval, Interval)
    assert (interval.chrom, interval.start, interval.end, interval.name, interval.score, interval.strand, interval.length) == ("chr4", 7, 11, "n", "3", "+", 4)
    path = tmp_path / "a.bed"
    bt = BedTool([interval]).saveas(path)
    assert str(BedTool(path)) == str(bt)
    assert bt.count() == 1
    assert str(BedTool("chr2\t1\t2\nchr1\t2\t3\n", from_string=True).sort()) == "chr1\t2\t3\nchr2\t1\t2\n"


def test_invalid_flag_combinations_are_rejected():
    with pytest.raises(ValueError):
        ours(A_TEXT).intersect(ours(B_TEXT), s=True, S=True)
    with pytest.raises(ValueError):
        ours(A_TEXT).intersect(ours(B_TEXT), f=0.0, r=True)
