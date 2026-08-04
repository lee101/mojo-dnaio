from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import ctypes

import pytest

import dnaio
from dnaio._lib import lib


def upstream_records(path):
    """Run the installed PyPI dnaio with this checkout removed from PYTHONPATH."""
    code = """import json, sys, dnaio
r=dnaio.FastqReader(sys.argv[1])
print(json.dumps({'two_headers': r.two_headers, 'records': [(x.name,x.sequence,x.qualities) for x in r]}))
"""
    env = os.environ.copy()
    env["PYTHONPATH"] = ""
    output = subprocess.check_output([sys.executable, "-c", code, str(path)], env=env, text=True)
    return json.loads(output)


@pytest.mark.parametrize("payload", [
    b"@one comment\nACGTN\n+\n!#%&'\n@two\nA\n+\nI\n",
    b"@one\r\nACGT\r\n+one\r\n!!!!\r\n@two\r\n\r\n+\r\n\r\n",
    b"@last\nAC\n+\n!!",  # dnaio deliberately accepts an absent final newline
])
def test_reader_matches_pypi_dnaio(tmp_path, payload):
    path = tmp_path / "reads.fastq"
    path.write_bytes(payload)
    expected = upstream_records(path)
    reader = dnaio.FastqReader(path)
    got = [(r.name, r.sequence, r.qualities) for r in reader]
    assert got == [tuple(x) for x in expected["records"]]
    assert reader.two_headers == expected["two_headers"]
    assert reader.number_of_records == len(expected["records"])


@pytest.mark.parametrize("payload, line", [
    (b"not-a-header\nA\n+\n!\n", 0),
    (b"@r\nAC\n+\n!\n", 3),
    (b"@r\nA\n+wrong\n!\n", 2),
    (b"@r\nA\n", 2),
])
def test_reader_rejects_invalid_fastq_at_matching_line(payload, line):
    with pytest.raises(dnaio.FastqFormatError) as error:
        dnaio.FastqReader(io.BytesIO(payload))
    assert error.value.line == line


def test_reader_rejects_non_ascii_fastq():
    with pytest.raises(dnaio.FastqFormatError) as error:
        dnaio.FastqReader(io.BytesIO(b"@r\nA\n+\n\xff\n"))
    assert error.value.line == 3


@pytest.mark.parametrize("two_headers", [False, True])
def test_writer_byte_parity_with_pypi_dnaio(tmp_path, two_headers):
    records = [dnaio.SequenceRecord("one note", "ACGT", "!#%&"), dnaio.SequenceRecord("two", "N", "I")]
    local = io.BytesIO()
    with dnaio.FastqWriter(local, two_headers=two_headers) as writer:
        for record in records:
            writer.write(record)
    expected_path = tmp_path / "upstream.fastq"
    code = """import sys, dnaio
with dnaio.FastqWriter(sys.argv[1], two_headers=sys.argv[2] == '1') as w:
    w.write(dnaio.SequenceRecord('one note', 'ACGT', '!#%&'))
    w.write(dnaio.SequenceRecord('two', 'N', 'I'))
"""
    env = os.environ.copy()
    env["PYTHONPATH"] = ""
    subprocess.check_call([sys.executable, "-c", code, str(expected_path), str(int(two_headers))], env=env)
    assert local.getvalue() == expected_path.read_bytes()


def test_mojo_formatter_simd_tail():
    name = b"n" * 151
    sequence = b"A" * 151
    qualities = b"I" * 151
    name_buffer = ctypes.create_string_buffer(name)
    sequence_buffer = ctypes.create_string_buffer(sequence)
    qualities_buffer = ctypes.create_string_buffer(qualities)
    target = ctypes.create_string_buffer(len(name) * 2 + len(sequence) + len(qualities) + 6)
    written = lib().md_format_fastq(ctypes.addressof(name_buffer), len(name),
                                    ctypes.addressof(sequence_buffer), len(sequence),
                                    ctypes.addressof(qualities_buffer), len(qualities), 1,
                                    ctypes.addressof(target), len(target))
    assert target.raw[:written] == b"@" + name + b"\n" + sequence + b"\n+" + name + b"\n" + qualities + b"\n"


def test_mojo_formatter_rejects_insufficient_destination():
    value = ctypes.create_string_buffer(b"A")
    target = ctypes.create_string_buffer(1)
    assert lib().md_format_fastq(ctypes.addressof(value), 1, ctypes.addressof(value), 1,
                                 ctypes.addressof(value), 1, 0, ctypes.addressof(target), 1) == -2


def test_mojo_scanner_simd_tail():
    name = b"n" * 151
    sequence = b"A" * 151
    payload = b"@" + name + b"\n" + sequence + b"\n+" + name + b"\n" + b"I" * 151
    reader = dnaio.FastqReader(io.BytesIO(payload))
    assert [(record.name, record.sequence, record.qualities) for record in reader] == [
        (name.decode(), sequence.decode(), "I" * 151)
    ]


def test_reader_preserves_custom_sequence_class():
    class CustomRecord:
        def __init__(self, name, sequence, qualities):
            self.values = name, sequence, qualities

    reader = dnaio.FastqReader(io.BytesIO(b"@r\nA\n+\n!\n"), sequence_class=CustomRecord)
    assert [record.values for record in reader] == [("r", "A", "!")]


def test_sequence_record_core_behavior_matches_upstream_surface():
    record = dnaio.SequenceRecord("read/1 comment", "ACGTRY", "abcdef")
    assert record.id == "read/1"
    assert record.comment == "comment"
    assert record.fastq_bytes(two_headers=True) == b"@read/1 comment\nACGTRY\n+read/1 comment\nabcdef\n"
    reversed_record = record.reverse_complement()
    assert (reversed_record.name, reversed_record.sequence, reversed_record.qualities) == ("read/1 comment", "RYACGT", "fedcba")
    assert dnaio.records_are_mates(dnaio.SequenceRecord("read/1", "A", "!"), dnaio.SequenceRecord("read/2", "A", "!"))


def test_open_api_roundtrip(tmp_path):
    path = tmp_path / "roundtrip.fastq"
    record = dnaio.SequenceRecord("r", "AC", "!!")
    with dnaio.open(path, mode="w", fileformat="fastq") as writer:
        writer.write(record)
    with dnaio.open(path, mode="r") as reader:
        assert list(reader) == [record]
    with dnaio.open(path, mode="a") as writer:
        writer.write(record)
    assert len(list(dnaio.open(path))) == 2
