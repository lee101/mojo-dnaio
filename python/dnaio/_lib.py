"""ctypes bridge for the Mojo FASTQ scanner."""

from __future__ import annotations

import ctypes
import os
import shutil
import subprocess

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LIB = os.environ.get("MOJO_DNAIO_LIB") or os.path.join(ROOT, "dist", "libmojo-dnaio.so")
I = ctypes.c_int64


def build() -> str:
    source = os.path.join(ROOT, "src", "fastq.mojo")
    if os.path.exists(LIB) and os.path.getmtime(LIB) >= os.path.getmtime(source):
        return LIB
    if os.environ.get("MOJO_DNAIO_LIB"):
        raise RuntimeError(f"MOJO_DNAIO_LIB does not exist or is stale: {LIB}")
    mojo = shutil.which("mojo")
    if not mojo:
        pixi = shutil.which("pixi") or os.path.expanduser("~/.pixi/bin/pixi")
        mojo = pixi
        command = [pixi, "run", "--manifest-path", os.path.join(ROOT, "pixi.toml"), "mojo"]
    else:
        command = [mojo]
    os.makedirs(os.path.dirname(LIB), exist_ok=True)
    proc = subprocess.run(command + ["build", "--emit", "shared-lib", source, "-o", LIB],
                          capture_output=True, text=True, timeout=1800)
    if proc.returncode or not os.path.exists(LIB):
        raise RuntimeError((proc.stderr or proc.stdout).strip()[:4000])
    return LIB


_cdll: ctypes.CDLL | None = None


def lib() -> ctypes.CDLL:
    global _cdll
    if _cdll is None:
        _cdll = ctypes.CDLL(build())
        _cdll.md_scan_fastq.argtypes = [I, I, I, I]
        _cdll.md_scan_fastq.restype = I
        _cdll.md_format_fastq.argtypes = [I] * 9
        _cdll.md_format_fastq.restype = I
    return _cdll


def bytes_array(value: bytes) -> np.ndarray:
    """A non-empty, C-contiguous byte array suitable for the C ABI."""
    return np.frombuffer(value if value else b"\0", dtype=np.uint8)


def format_fastq(name: str, sequence: str, qualities: str, two_headers: bool) -> np.ndarray:
    """Format one validated record into a caller-owned uint8 array.

    The local NumPy arrays remain strongly referenced for the entire foreign call;
    their dtypes and contiguity are fixed here rather than assumed by the C ABI.
    """
    name_bytes = bytes_array(name.encode("ascii"))
    sequence_bytes = bytes_array(sequence.encode("ascii"))
    quality_bytes = bytes_array(qualities.encode("ascii"))
    size = len(name_bytes) * (2 if two_headers else 1) + len(sequence_bytes) + len(quality_bytes) + 6
    target = np.empty(size, dtype=np.uint8)
    written = lib().md_format_fastq(
        name_bytes.ctypes.data, len(name), sequence_bytes.ctypes.data, len(sequence),
        quality_bytes.ctypes.data, len(qualities), int(two_headers), target.ctypes.data, target.size,
    )
    if written < 0 or written != size:
        raise RuntimeError("Mojo FASTQ formatter rejected validated input")
    return target
