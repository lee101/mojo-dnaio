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

`src/fastq.mojo` is one compilation unit exported as `dist/libmojo-dnaio.so`. The ctypes wrapper supplies a non-null, zero-copy NumPy view of the input and contiguous `int64` offset storage, keeping both buffers alive for the complete foreign call. Mojo uses one SIMD pass to find line endings and validate ASCII, then writes six offsets per record (name, sequence, qualities). Python decodes only those slices when iteration materializes slotted `SequenceRecord`s. The normal writer creates one final `bytes` object per record; the exported SIMD Mojo formatter remains available for callers that already own byte buffers.

## Benchmarks

Measured by `pixi run bench` on this machine on 2026-08-24, with a 101.7 MB synthetic FASTQ containing 300,000 records of 150 bases. The reported time is the best of three iterations.

| operation | mojo-dnaio | dnaio 1.2.4 | speedup |
| --- | ---: | ---: | ---: |
| parse and materialize records | 784.16 ms | 228.60 ms | 0.29× |
| format and write records | 239.44 ms | 135.41 ms | 0.57× |

These end-to-end operations remain slower than upstream `dnaio` 1.2.4. Its Cython extension constructs Python unicode objects in-process; this standalone port materializes them in Python. The scan is sequential because record boundaries are variable, and formatting writes into one ordered stream, so thread launch and synchronization would add overhead rather than expose useful independent work. No GPU path is included: scanning and formatting are streaming byte operations with far less than two arithmetic operations per byte moved, so device transfer would dominate.

## Validation

`pixi run test` asserts parser records, repeated-header detection, CRLF behavior, absent-final-newline behavior, error locations, and writer bytes against the PyPI `dnaio` 1.2.4 package installed in the Pixi environment.
