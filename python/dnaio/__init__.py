"""Mojo-accelerated compatible subset of dnaio: single-end FASTQ I/O."""

from __future__ import annotations

import os

from ._core import SequenceRecord, record_names_match, records_are_mates
from .exceptions import FastqFormatError, FileFormatError, UnknownFileFormat
from .readers import FastqReader
from .writers import FastqWriter

Sequence = SequenceRecord
__version__ = "0.1.0"


def open(*files, file1=None, file2=None, fileformat=None, interleaved=False, mode="r",
         qualities=None, opener=open, compression_level=1, open_threads=0, **_kwargs):
    """Open an uncompressed, single-end FASTQ source for reading, writing, or appending."""
    if file1 is not None or file2 is not None:
        if files: raise ValueError("file1 and file2 cannot be combined with positional files")
        files = tuple(x for x in (file1, file2) if x is not None)
    if len(files) != 1 or interleaved:
        raise NotImplementedError("mojo-dnaio currently supports one single-end FASTQ file")
    if fileformat not in (None, "fastq", "fq"):
        raise UnknownFileFormat("Only FASTQ is covered by mojo-dnaio")
    if mode == "r": return FastqReader(files[0], opener=opener)
    if mode in ("w", "a"):
        if isinstance(files[0], (str, bytes, os.PathLike)):
            file_mode = "ab" if mode == "a" else "wb"
            return FastqWriter(opener(files[0], file_mode), _close_file=True)
        return FastqWriter(files[0])
    raise ValueError("Mode must be 'r', 'w' or 'a'")
