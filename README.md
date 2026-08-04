# mojo-pybedtools

`mojo-pybedtools` is a standalone, in-memory port of the compute-heavy BED interval
arithmetic at the heart of [pybedtools](https://daler.github.io/pybedtools/).  It keeps
the familiar `BedTool` and `Interval` surface for the covered operations, but runs the
sweep and merge kernels in Mojo instead of creating temporary files and invoking the
`bedtools` executable.

## Covered subset

`BedTool.intersect()` supports the default clipped result and `wa`, `wb`, `u`, `v`,
`c`, `wo`, `f`, `F`, `r`, `e`, `s`, and `S`. `BedTool.subtract()` supports `A` and the
same overlap filters; `BedTool.coverage()`, `BedTool.merge(d=...)`, `sort()`, `count()`,
`saveas()`, `Interval`, and `create_interval_from_list()` are also implemented.

This is intentionally not the full pybedtools wrapper. BAM/VCF/GFF parsing, temporary
file workflows, streaming, genome-file ordering, `closest`, `window`, `slop`,
`shuffle`, and the many less-common bedtools flags are not covered. Unsupported options
raise `NotImplementedError` rather than silently producing a different genomic result.

## Install and use

```bash
pixi install
pixi run build
pixi run test
```

```python
from mojo_pybedtools import BedTool

genes = BedTool("""chr1\t10\t40\tgeneA\t0\t+
chr1\t80\t100\tgeneB\t0\t-
""", from_string=True)
peaks = BedTool("""chr1\t20\t30\tpeak1\t0\t+
chr1\t95\t110\tpeak2\t0\t-
""", from_string=True)

print(genes.intersect(peaks, wa=True, wb=True), end="")
# chr1    10    40    geneA    0    +    chr1    20    30    peak1    0    +
# chr1    80    100   geneB    0    -    chr1    95    110   peak2    0    -

print(genes.coverage(peaks), end="")
# ... count, covered bases, A length, and covered fraction are appended per record.
```

## Benchmarks

Measured with `pixi run bench` on this machine (`x86_64`), Python 3.13.14. Each input
has 25,000 random BED intervals across 22
chromosomes; timings include materialising the result text and are the best of three.

| kernel | mojo-pybedtools | pybedtools/bedtools | speedup |
| --- | ---: | ---: | ---: |
| intersect -c | 88.64 ms | 232.37 ms | 2.62x |
| coverage | 170.67 ms | 366.37 ms | 2.15x |
| merge | 97.24 ms | 249.23 ms | 2.56x |

These are real local measurements, not projections. The Python-facing parsing and BED
record construction remain visible at this size. The interval sweeps are branch-heavy,
memory-bound operations, so no GPU path is provided: host/device transfer costs exceed
their arithmetic work.

## How it works

The wrapper normalises BED records into contiguous `int64` columns for chromosome code,
start, end, and strand. It caches columns sorted for each input pairing, passes their raw
addresses to a single Mojo shared library through `ctypes`, and maps output indices back
to the original BED fields. The Mojo kernels use a sweep frontier over sorted B starts:
expired intervals are skipped, candidate overlaps are scanned only until A's end, and
all coordinates use BED's half-open `[start, end)` convention. No compiled code owns
Python memory, so there is no cross-language allocation or lifetime management.

`tests/test_parity.py` compares all covered operations and flags directly with
pybedtools 0.12.0 plus bedtools 2.31.1 installed by Pixi.

MIT. See [LICENSE](LICENSE).
