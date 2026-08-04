from __future__ import annotations

import os

from ._core import SequenceRecord
from ._lib import format_fastq


class FastqWriter:
    """Write FASTQ records; Mojo performs the byte-level record formatting."""
    def __init__(self, file, *, two_headers=False, opener=open, _close_file=None):
        self._two_headers = two_headers
        self._close_on_exit = False
        if isinstance(file, (str, os.PathLike)):
            self._file = opener(file, "wb")
            self._close_on_exit = True
        else:
            self._file = file
            self._close_on_exit = bool(_close_file)

    def write(self, record: SequenceRecord) -> None:
        if not isinstance(record, SequenceRecord):
            raise TypeError("record must be a SequenceRecord")
        if record.qualities is None:
            raise ValueError("Cannot write FASTQ record without qualities")
        self._file.write(format_fastq(record.name, record.sequence, record.qualities, self._two_headers))

    def writeseq(self, name: str, sequence: str, qualities: str) -> None:
        self.write(SequenceRecord(name, sequence, qualities))

    def close(self):
        if self._close_on_exit: self._file.close()

    def __enter__(self): return self
    def __exit__(self, *args): self.close()
