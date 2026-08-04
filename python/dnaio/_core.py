"""Python objects matching dnaio's FASTQ-facing record API."""

from __future__ import annotations

from dataclasses import dataclass


def _ascii(field: str, value: str) -> str:
    if type(value) is not str:
        raise TypeError(f"{field} should be of type str, got {type(value)}")
    try:
        value.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError(f"'{field}' in sequence file must be ASCII encoded, but found '{value[exc.start:exc.end]}' at index {exc.start}") from None
    return value


@dataclass(eq=True)
class SequenceRecord:
    name: str
    sequence: str
    qualities: str | None = None

    def __post_init__(self):
        self.name = _ascii("name", self.name)
        self.sequence = _ascii("sequence", self.sequence)
        if self.qualities is not None:
            self.qualities = _ascii("qualities", self.qualities)
            if len(self.sequence) != len(self.qualities):
                raise ValueError(f"In read named {self.name!r}: length of quality sequence ({len(self.qualities)}) and length of read ({len(self.sequence)}) do not match")

    @classmethod
    def _from_validated_fastq(cls, name: str, sequence: str, qualities: str) -> "SequenceRecord":
        record = cls.__new__(cls)
        record.name = name
        record.sequence = sequence
        record.qualities = qualities
        return record

    @property
    def id(self) -> str:
        return self.name.split(None, 1)[0]

    @property
    def comment(self) -> str | None:
        parts = self.name.split(None, 1)
        return parts[1] if len(parts) == 2 else None

    def fastq_bytes(self, two_headers: bool = False) -> bytes:
        if self.qualities is None:
            raise ValueError("Cannot write FASTQ record without qualities")
        plus = self.name if two_headers else ""
        return f"@{self.name}\n{self.sequence}\n+{plus}\n{self.qualities}\n".encode("ascii")

    def reverse_complement(self) -> "SequenceRecord":
        table = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
        return SequenceRecord(self.name, self.sequence.translate(table)[::-1],
                              None if self.qualities is None else self.qualities[::-1])

    def is_mate(self, other: "SequenceRecord") -> bool:
        return records_are_mates(self, other)


def records_are_mates(record1: SequenceRecord, record2: SequenceRecord) -> bool:
    left, right = record1.id, record2.id
    if left == right:
        return True
    return len(left) == len(right) and left[:-1] == right[:-1] and left[-1:] in "123" and right[-1:] in "123" and left[-1] != right[-1]


record_names_match = records_are_mates
