from __future__ import annotations

import os
from typing import BinaryIO, Iterator

import numpy as np

from ._core import SequenceRecord
from ._lib import bytes_array, lib
from .exceptions import FastqFormatError


class FastqReader:
    """Read single-line FASTQ records using Mojo to scan and validate bytes."""
    delivers_qualities = True

    def __init__(self, file, *, sequence_class=SequenceRecord, buffer_size=128 * 1024,
                 opener=open, _close_file=None):
        self.sequence_class = sequence_class
        self.buffer_size = buffer_size
        self._close_on_exit = False
        if isinstance(file, (str, os.PathLike)):
            self._file = opener(file, "rb")
            self._close_on_exit = True
        else:
            self._file = file
            self._close_on_exit = bool(_close_file)
        raw = self._file.read()
        if not isinstance(raw, bytes):
            raise TypeError("self.file is not a binary file reader.")
        self._data = raw
        self._offsets = self._scan(raw)
        self.number_of_records = len(self._offsets)
        self.two_headers = bool(self.number_of_records and self._has_two_headers())

    def _scan(self, raw: bytes):
        if not raw:
            return []
        data = bytes_array(raw)
        # Every complete FASTQ record has four line endings. `bytes.count` is a
        # fast C-level upper bound and avoids reserving six int64s per input byte.
        capacity = raw.count(b"\n") // 4 + 1
        fields = np.empty(capacity * 6, dtype=np.int64)
        status = lib().md_scan_fastq(data.ctypes.data, len(raw), fields.ctypes.data, capacity)
        if status < 0:
            line = -status - 1
            messages = {0: "Line expected to start with '@'", 1: "Premature end of file encountered.",
                        2: "Line expected to start with '+'", 3: "Length of sequence and qualities differ"}
            raise FastqFormatError(messages.get(line % 4, "Invalid FASTQ record"), line)
        return fields[:status * 6].reshape(status, 6).tolist()

    def _has_two_headers(self) -> bool:
        nstart, nlen, sstart, slen, qstart, qlen = self._offsets[0]
        plus_start = sstart + slen
        # Account for CRLF on the sequence line.
        if self._data[plus_start:plus_start + 2] == b"\r\n": plus_start += 2
        else: plus_start += 1
        return self._data[plus_start + 1: qstart].rstrip(b"\r\n") != b""

    def __iter__(self) -> Iterator[SequenceRecord]:
        if self.sequence_class is SequenceRecord:
            for nstart, nlen, sstart, slen, qstart, qlen in self._offsets:
                yield SequenceRecord._from_validated_fastq(
                    self._data[nstart:nstart+nlen].decode("ascii"),
                    self._data[sstart:sstart+slen].decode("ascii"),
                    self._data[qstart:qstart+qlen].decode("ascii"))
            return
        for nstart, nlen, sstart, slen, qstart, qlen in self._offsets:
            yield self.sequence_class(self._data[nstart:nstart+nlen].decode("ascii"),
                                      self._data[sstart:sstart+slen].decode("ascii"),
                                      self._data[qstart:qstart+qlen].decode("ascii"))

    def close(self):
        if self._close_on_exit and self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self):
        if self._file is None: raise ValueError("I/O operation on closed BinaryFileReader")
        return self

    def __exit__(self, *args): self.close()
