# mojo-dnaio

A standalone Mojo port of the compute-heavy core of [dnaio](https://github.com/marcelm/dnaio): validating and scanning single-line FASTQ records. It keeps the covered Python names and signatures familiar, while returning normal Python `SequenceRecord` objects.

Covered: `SequenceRecord`, `FastqReader`, `FastqWriter`, `records_are_mates`, and `dnaio.open()` for one uncompressed single-end FASTQ stream. The reader accepts LF and CRLF files and an absent final newline, validates headers, repeated `+name` headers, and sequence/quality lengths. Not covered: FASTA, gzip/xopen, BAM, paired/interleaved/multiple files, and chunk APIs.

## Install and use

```bash
pixi install
pixi run build
pixi run python - <<'PY'
import dnaio

with dnaio.open("reads.fastq", mode="w", fileformat="fastq") as writer:
    writer.write(dnaio.SequenceRecord("read-1", "ACGT", "!!!!"))
with dnaio.open("reads.fastq") as reader:
    print(next(iter(reader)))
PY
```

## How it works

`src/fastq.mojo` is one compilation unit exported as `dist/libmojo-dnaio.so`. The ctypes wrapper supplies non-null, C-contiguous NumPy `uint8` input/output buffers and `int64` offset storage, keeping each buffer alive for the complete foreign call. Mojo scans the input once, validates each four-line record, and writes six offsets per record (name, sequence, qualities). Python decodes only those slices when iteration materializes `SequenceRecord`s. The writer uses the reciprocal Mojo formatter into a capacity-checked caller-owned `uint8` buffer.

## Benchmarks

Measured by `pixi run bench` on this machine on 2026-08-04, with a 101.7 MB synthetic FASTQ containing 300,000 records of 150 bases. The reported time is the best of three iterations.

| operation | mojo-dnaio | dnaio 1.2.4 | speedup |
| --- | ---: | ---: | ---: |
| parse and materialize records | 892.91 ms | 269.93 ms | 0.30× |
| format and write records | 4065.25 ms | 137.86 ms | 0.03× |

These end-to-end operations are slower than upstream `dnaio` 1.2.4. Its Cython extension constructs Python unicode objects in-process; this standalone port materializes them in Python. Formatting is particularly slower because the safe ctypes bridge allocates a caller-owned NumPy buffer per record. SIMD accelerates the Mojo byte loops, but the benchmark includes the Python-side work. No GPU path is included: scanning and formatting are streaming byte operations with too little arithmetic intensity to overcome device-transfer overhead.

## Validation

`pixi run test` asserts parser records, repeated-header detection, CRLF behavior, absent-final-newline behavior, error locations, and writer bytes against the PyPI `dnaio` 1.2.4 package installed in the Pixi environment.
