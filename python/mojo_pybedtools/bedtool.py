"""A compact, in-memory subset of :mod:`pybedtools` powered by Mojo."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np

from ._lib import lib


@dataclass(frozen=True, slots=True)
class Interval:
    """A BED record with pybedtools-compatible core attributes."""

    fields: tuple[str, ...]

    def __init__(self, *fields: object):
        if len(fields) == 1 and isinstance(fields[0], str):
            fields = tuple(fields[0].rstrip("\n").split("\t"))
        if len(fields) < 3:
            raise ValueError("a BED interval needs chrom, start, and end")
        chrom, start, end = str(fields[0]), int(fields[1]), int(fields[2])
        if start < 0 or end < start:
            raise ValueError("BED coordinates must satisfy 0 <= start <= end")
        object.__setattr__(self, "fields", (chrom, str(start), str(end)) + tuple(map(str, fields[3:])))

    @classmethod
    def _from_fields(cls, fields: tuple[str, ...]) -> "Interval":
        interval = object.__new__(cls)
        object.__setattr__(interval, "fields", fields)
        return interval

    @property
    def chrom(self) -> str:
        return self.fields[0]

    @property
    def start(self) -> int:
        return int(self.fields[1])

    @property
    def end(self) -> int:
        return int(self.fields[2])

    @property
    def name(self) -> str:
        return self.fields[3] if len(self.fields) > 3 else "."

    @property
    def score(self) -> str:
        return self.fields[4] if len(self.fields) > 4 else "."

    @property
    def strand(self) -> str:
        return self.fields[5] if len(self.fields) > 5 else "."

    @property
    def length(self) -> int:
        return self.end - self.start

    def __str__(self) -> str:
        return "\t".join(self.fields) + "\n"


def _coerce_item(item: object) -> Interval:
    if isinstance(item, Interval):
        return item
    if isinstance(item, str):
        return Interval(item)
    return Interval(*item)  # type: ignore[arg-type]


class BedTool:
    """In-memory BED collection for the interval arithmetic subset."""

    def __init__(self, fn: str | Path | Iterable[object] | "BedTool", from_string: bool = False, **_: object):
        if isinstance(fn, BedTool):
            self._records = fn._records
        elif isinstance(fn, (str, Path)):
            text = str(fn) if from_string else Path(fn).read_text() if Path(fn).exists() else str(fn)
            self._records = tuple(Interval(line) for line in text.splitlines() if line and not line.startswith(("#", "track", "browser")))
        else:
            self._records = tuple(_coerce_item(x) for x in fn)
        self._chrom = np.asarray([record.chrom for record in self._records])
        self._start = np.fromiter((record.start for record in self._records), dtype=np.int64, count=len(self._records))
        self._end = np.fromiter((record.end for record in self._records), dtype=np.int64, count=len(self._records))
        self._strand = np.fromiter(({"+": 1, "-": -1}.get(record.strand, 0) for record in self._records), dtype=np.int64, count=len(self._records))
        self._encoded_cache: dict[BedTool | None, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}

    @classmethod
    def _from_records(
        cls,
        records: tuple[Interval, ...],
        chrom: np.ndarray,
        start: np.ndarray,
        end: np.ndarray,
        strand: np.ndarray,
    ) -> "BedTool":
        tool = object.__new__(cls)
        tool._records = records
        tool._chrom = chrom
        tool._start = start
        tool._end = end
        tool._strand = strand
        tool._encoded_cache = {}
        return tool

    def __iter__(self) -> Iterator[Interval]:
        return iter(self._records)

    def __len__(self) -> int:
        return len(self._records)

    def __str__(self) -> str:
        return "".join("\t".join(record.fields) + "\n" for record in self._records)

    def __repr__(self) -> str:
        return f"BedTool({len(self)} intervals)"

    def count(self) -> int:
        return len(self)

    def saveas(self, fn: str | Path | None = None, **_: object) -> "BedTool":
        if fn is not None:
            Path(fn).write_text(str(self))
        return self

    def sort(self, **_: object) -> "BedTool":
        return BedTool(sorted(self._records, key=lambda x: (x.chrom, x.start, x.end, x.fields)))

    def _encoded(self, other: "BedTool | None" = None):
        cached = self._encoded_cache.get(other)
        if cached is not None:
            return cached, None if other is None else other._encoded_cache[self]
        chrom = self._chrom if other is None else np.concatenate((self._chrom, other._chrom))
        chroms = np.unique(chrom)

        def encode(tool: "BedTool"):
            codes = np.searchsorted(chroms, tool._chrom)
            order = np.lexsort((np.arange(len(tool), dtype=np.int64), tool._end, tool._start, codes))
            return (
                order,
                np.ascontiguousarray(codes[order], dtype=np.int64),
                np.ascontiguousarray(tool._start[order]),
                np.ascontiguousarray(tool._end[order]),
                np.ascontiguousarray(tool._strand[order]),
            )

        encoded = encode(self)
        self._encoded_cache[other] = encoded
        if other is None:
            return encoded, None
        other_encoded = encode(other)
        other._encoded_cache[self] = other_encoded
        return encoded, other_encoded

    def _counts(self, b: "BedTool", f: float, F: float, r: bool, e: bool, s: bool, S: bool) -> np.ndarray:
        (ao, ac, ast, aen, ass), (bo, bc, bst, ben, bss) = self._encoded(b)  # type: ignore[misc]
        counts = np.zeros(len(self), dtype=np.int64)
        if len(ao) and len(bo):
            sorted_counts = np.empty(len(ao), dtype=np.int64)
            args = [x.ctypes.data for x in (ac, ast, aen, ass, bc, bst, ben, bss)] + [len(ao), len(bo), float(f), float(F), int(r), int(e), 1 if s else 2 if S else 0, sorted_counts.ctypes.data]
            lib().mpbt_intersect_counts(*args)
            counts[ao] = sorted_counts
        return counts

    def _pairs(self, b: "BedTool", f: float = 0.0, F: float = 0.0, r: bool = False, e: bool = False, s: bool = False, S: bool = False):
        if s and S:
            raise ValueError("-s and -S are mutually exclusive")
        if r and not f:
            raise ValueError("-r requires f")
        (ao, ac, ast, aen, ass), (bo, bc, bst, ben, bss) = self._encoded(b)  # type: ignore[misc]
        if not len(ao) or not len(bo):
            return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
        args = [x.ctypes.data for x in (ac, ast, aen, ass, bc, bst, ben, bss)] + [len(ao), len(bo), float(f), float(F), int(r), int(e), 1 if s else 2 if S else 0]
        n = lib().mpbt_intersect_count(*args)
        left, right = np.empty(n, dtype=np.int64), np.empty(n, dtype=np.int64)
        written = lib().mpbt_intersect_pairs(*args, left.ctypes.data, right.ctypes.data)
        assert written == n
        return ao[left], bo[right]

    def intersect(self, b: "BedTool | str | Path | Iterable[object]", wa: bool = False, wb: bool = False, u: bool = False, c: bool = False, v: bool = False, f: float = 0.0, F: float = 0.0, r: bool = False, e: bool = False, s: bool = False, S: bool = False, wo: bool = False, **kwargs: object) -> "BedTool":
        """Return overlaps using the common ``bedtools intersect`` flags."""
        if kwargs:
            raise NotImplementedError(f"unsupported intersect options: {', '.join(kwargs)}")
        if sum(bool(x) for x in (u, c, v)) > 1:
            raise ValueError("u, c, and v are mutually exclusive")
        other = b if isinstance(b, BedTool) else BedTool(b)
        if c:
            counts = self._counts(other, f, F, r, e, s, S)
            records = tuple(
                Interval._from_fields(a.fields + (str(int(counts[i])),))
                for i, a in enumerate(self._records)
            )
            return BedTool._from_records(records, self._chrom, self._start, self._end, self._strand)
        left, right = self._pairs(other, f, F, r, e, s, S)
        hit = np.zeros(len(self), dtype=bool)
        hit[left] = True
        if u:
            return BedTool([a for i, a in enumerate(self._records) if hit[i]])
        if v:
            return BedTool([a for i, a in enumerate(self._records) if not hit[i]])
        pairs = sorted(zip(left.tolist(), right.tolist()))
        result: list[Interval] = []
        for i, j in pairs:
            a, bb = self._records[i], other._records[j]
            if wo:
                result.append(Interval(*a.fields, *bb.fields, min(a.end, bb.end) - max(a.start, bb.start)))
            elif wa and wb:
                result.append(Interval(*a.fields, *bb.fields))
            elif wa:
                result.append(a)
            elif wb:
                result.append(Interval(a.chrom, max(a.start, bb.start), min(a.end, bb.end), *a.fields[3:], *bb.fields))
            else:
                result.append(Interval(a.chrom, max(a.start, bb.start), min(a.end, bb.end), *a.fields[3:]))
        return BedTool(result)

    def merge(self, d: int = 0, **kwargs: object) -> "BedTool":
        """Merge overlapping/book-ended intervals, as ``bedtools merge -d``."""
        if kwargs:
            raise NotImplementedError(f"unsupported merge options: {', '.join(kwargs)}")
        (order, chrom, start, end, _), _ = self._encoded()
        # Mojo pointers are non-nullable.  NumPy does not promise a meaningful data
        # address for an empty array, so keep the empty case on the Python side.
        if not len(order):
            return BedTool(())
        dst_c = np.empty(len(order), dtype=np.int64)
        dst_s = np.empty(len(order), dtype=np.int64)
        dst_e = np.empty(len(order), dtype=np.int64)
        n = lib().mpbt_merge(chrom.ctypes.data, start.ctypes.data, end.ctypes.data, len(order), int(d), dst_c.ctypes.data, dst_s.ctypes.data, dst_e.ctypes.data)
        labels = np.unique(self._chrom)
        merged_chrom = labels[dst_c[:n]]
        merged_start = dst_s[:n]
        merged_end = dst_e[:n]
        records = tuple(
            Interval._from_fields((chrom_name, str(start_value), str(end_value)))
            for chrom_name, start_value, end_value in zip(
                merged_chrom.tolist(), merged_start.tolist(), merged_end.tolist()
            )
        )
        return BedTool._from_records(
            records,
            merged_chrom,
            merged_start,
            merged_end,
            np.zeros(n, dtype=np.int64),
        )

    def subtract(self, b: "BedTool | str | Path | Iterable[object]", A: bool = False, f: float = 0.0, F: float = 0.0, r: bool = False, e: bool = False, s: bool = False, S: bool = False, **kwargs: object) -> "BedTool":
        """Remove overlapping pieces of A; ``A=True`` drops any hit record."""
        if kwargs:
            raise NotImplementedError(f"unsupported subtract options: {', '.join(kwargs)}")
        other = b if isinstance(b, BedTool) else BedTool(b)
        left, right = self._pairs(other, f, F, r, e, s, S)
        if A:
            hit = np.zeros(len(self), dtype=bool)
            hit[left] = True
            return BedTool([a for i, a in enumerate(self._records) if not hit[i]])
        by_a: dict[int, list[Interval]] = {}
        for i, j in zip(left.tolist(), right.tolist()):
            by_a.setdefault(i, []).append(other._records[j])
        result: list[Interval] = []
        for i, a in enumerate(self._records):
            pos = a.start
            for bb in sorted(by_a.get(i, ()), key=lambda x: (x.start, x.end)):
                if bb.start > pos:
                    result.append(Interval(a.chrom, pos, min(bb.start, a.end), *a.fields[3:]))
                pos = max(pos, bb.end)
                if pos >= a.end:
                    break
            if pos < a.end:
                result.append(Interval(a.chrom, pos, a.end, *a.fields[3:]))
        return BedTool(result)

    def coverage(self, b: "BedTool | str | Path | Iterable[object]", **kwargs: object) -> "BedTool":
        """Append count, covered bases, length, and fraction as ``bedtools coverage``."""
        if kwargs:
            raise NotImplementedError(f"unsupported coverage options: {', '.join(kwargs)}")
        other = b if isinstance(b, BedTool) else BedTool(b)
        (ao, ac, ast, aen, _), (_, bc, bst, ben, _) = self._encoded(other)  # type: ignore[misc]
        bases, counts = np.zeros(len(ao), dtype=np.int64), np.zeros(len(ao), dtype=np.int64)
        if len(ao) and len(bc):
            lib().mpbt_coverage(ac.ctypes.data, ast.ctypes.data, aen.ctypes.data, bc.ctypes.data, bst.ctypes.data, ben.ctypes.data, ao.ctypes.data, len(ao), len(bc), bases.ctypes.data, counts.ctypes.data)
        lengths = (self._end - self._start).tolist()
        records = tuple(
            Interval._from_fields(
                a.fields + (
                    str(count),
                    str(base),
                    str(length),
                    f"{base / length if length else 0.0:.7f}",
                )
            )
            for a, count, base, length in zip(
                self._records, counts.tolist(), bases.tolist(), lengths
            )
        )
        return BedTool._from_records(records, self._chrom, self._start, self._end, self._strand)


def create_interval_from_list(fields: Sequence[object]) -> Interval:
    return Interval(*fields)


def example_bedtool(filename: str) -> BedTool:
    return BedTool(filename)
