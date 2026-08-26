"""ctypes binding for the standalone Mojo interval kernels."""

from __future__ import annotations

import ctypes
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LIB = os.path.join(ROOT, "dist", "libmojo-pybedtools.so")
I = ctypes.c_int64
F = ctypes.c_double

_SIGNATURES = {
    "mpbt_intersect_count": ([I] * 10 + [F, F] + [I] * 3, I),
    "mpbt_intersect_pairs": ([I] * 10 + [F, F] + [I] * 5, I),
    "mpbt_intersect_counts": ([I] * 10 + [F, F] + [I] * 4, None),
    "mpbt_coverage": ([I] * 11, None),
    "mpbt_merge": ([I] * 8, I),
}
_cached: ctypes.CDLL | None = None


def build() -> str:
    source = os.path.join(ROOT, "src", "capi.mojo")
    if not os.path.exists(LIB) or os.path.getmtime(LIB) < os.path.getmtime(source):
        subprocess.run(["bash", os.path.join(ROOT, "build", "build.sh")], check=True, timeout=1800)
    return LIB


def lib() -> ctypes.CDLL:
    global _cached
    if _cached is None:
        _cached = ctypes.CDLL(build())
        for name, (args, result) in _SIGNATURES.items():
            fn = getattr(_cached, name)
            fn.argtypes, fn.restype = args, result
    return _cached
