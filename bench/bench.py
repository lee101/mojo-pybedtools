"""Measured Mojo-versus-bedtools interval arithmetic benchmark."""

from __future__ import annotations

import os
import platform
import sys
import time

import numpy as np
import pybedtools

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
from mojo_pybedtools import BedTool  # noqa: E402


def render(chrom: np.ndarray, start: np.ndarray, end: np.ndarray) -> str:
    return "".join(f"chr{int(c)}\t{int(s)}\t{int(e)}\n" for c, s, e in zip(chrom, start, end))


def best(fn, reps: int = 3) -> float:
    result = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        result = min(result, time.perf_counter() - t0)
    return result


def row(kernel: str, mojo: float, reference: float) -> None:
    print(f"| {kernel} | {mojo * 1e3:.2f} ms | {reference * 1e3:.2f} ms | {reference / mojo:.2f}x |")


def main() -> None:
    rng = np.random.default_rng(7)
    n = 25_000
    chrom_a = rng.integers(1, 23, n)
    start_a = rng.integers(0, 20_000_000, n)
    end_a = start_a + rng.integers(100, 1_000, n)
    chrom_b = rng.integers(1, 23, n)
    start_b = rng.integers(0, 20_000_000, n)
    end_b = start_b + rng.integers(100, 1_000, n)
    a_text, b_text = render(chrom_a, start_a, end_a), render(chrom_b, start_b, end_b)
    a, b = BedTool(a_text, from_string=True), BedTool(b_text, from_string=True)
    ua, ub = pybedtools.BedTool(a_text, from_string=True), pybedtools.BedTool(b_text, from_string=True)
    # Build and load before timing.  Each operation is materialised, matching bedtools.
    a.intersect(b, c=True)
    print(f"Machine: {platform.processor() or platform.machine()} | Python {platform.python_version()} | n={n:,} per input")
    print("| kernel | mojo-pybedtools | pybedtools/bedtools | speedup |")
    print("| --- | ---: | ---: | ---: |")
    row("intersect -c", best(lambda: str(a.intersect(b, c=True))), best(lambda: str(ua.intersect(ub, c=True))))
    row("coverage", best(lambda: str(a.coverage(b))), best(lambda: str(ua.coverage(ub))))
    row("merge", best(lambda: str(a.merge())), best(lambda: str(ua.sort().merge())))


if __name__ == "__main__":
    main()
