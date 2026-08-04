"""Benchmark the Mojo scanner and formatter against installed dnaio."""

from __future__ import annotations

import io
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))

import dnaio


def timeit(fn, reps=3):
    best = float("inf")
    result = None
    for _ in range(reps):
        start = time.perf_counter()
        result = fn()
        best = min(best, time.perf_counter() - start)
    return best, result


def make_fastq(records=300_000, length=150):
    sequence = "ACGT" * (length // 4) + "A" * (length % 4)
    qualities = "I" * length
    line = f"@instrument:1:flow:1:1:1:1 comment\n{sequence}\n+\n{qualities}\n".encode()
    return line * records


def row(name, mojo, upstream):
    print(f"| {name} | {mojo * 1e3:.2f} ms | {upstream * 1e3:.2f} ms | {upstream / mojo:.2f}x |")


def main():
    payload = make_fastq()
    print(f"payload: {len(payload) / 1e6:.1f} MB, 300,000 records x 150 bases")
    print("| operation | mojo-dnaio | dnaio 1.2.4 | speedup |")
    print("| --- | ---: | ---: | ---: |")
    local, records = timeit(lambda: list(dnaio.FastqReader(io.BytesIO(payload))))

    import subprocess
    env = os.environ.copy(); env["PYTHONPATH"] = ""
    def upstream_time(body):
        code = f"""import io, sys, time, dnaio
data = sys.stdin.buffer.read()
records = list(dnaio.FastqReader(io.BytesIO(data)))
best = float('inf')
for _ in range(3):
    start = time.perf_counter()
    {body}
    best = min(best, time.perf_counter() - start)
print(best)
"""
        return float(subprocess.check_output([sys.executable, "-c", code], input=payload, env=env).decode())
    upstream = upstream_time("parsed = list(dnaio.FastqReader(io.BytesIO(data)))")
    row("parse and materialize records", local, upstream)

    def local_write():
        stream = io.BytesIO()
        with dnaio.FastqWriter(stream) as writer:
            for record in records:
                writer.write(record)
        return stream.getvalue()
    local_write_s, _ = timeit(local_write)
    upstream_write_s = upstream_time("stream = io.BytesIO(); [stream.write(x.fastq_bytes()) for x in records]")
    row("format and write records", local_write_s, upstream_write_s)


if __name__ == "__main__":
    main()
