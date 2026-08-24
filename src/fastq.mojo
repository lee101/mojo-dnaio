"""FASTQ record scanner and formatter used through the small C ABI below."""

from std.sys import simd_width_of

comptime BPtr = UnsafePointer[UInt8, AnyOrigin[mut=True]]
comptime IPtr = UnsafePointer[Int, AnyOrigin[mut=True]]


def ascii_line_end(data: BPtr, n: Int, start: Int) -> Int:
    """Return the LF offset, -1 at EOF, or -2 for non-ASCII input."""
    comptime W = simd_width_of[DType.uint8]()
    var i = start
    while i + W <= n:
        var values = data.load[width=W](i)
        if values.eq(SIMD[DType.uint8, W](UInt8(10))).reduce_or()[0]:
            for j in range(W):
                if data[i + j] == UInt8(10):
                    return i + j
                if data[i + j] > UInt8(127):
                    return -2
        elif not values.le(SIMD[DType.uint8, W](UInt8(127))).reduce_and()[0]:
            return -2
        i += W
    while i < n:
        if data[i] == UInt8(10):
            return i
        if data[i] > UInt8(127):
            return -2
        i += 1
    return -1


def trimmed_length(data: BPtr, start: Int, end: Int) -> Int:
    if end > start and data[end - 1] == UInt8(13):
        return end - start - 1
    return end - start


def equal_bytes(data: BPtr, left: Int, right: Int, width: Int) -> Bool:
    comptime W = simd_width_of[DType.uint8]()
    var i = 0
    while i + W <= width:
        if data.load[width=W](left + i).ne(data.load[width=W](right + i)).reduce_or()[0]:
            return False
        i += W
    while i < width:
        if data[left + i] != data[right + i]:
            return False
        i += 1
    return True


def scan_fastq(data: BPtr, n: Int, fields: IPtr, capacity: Int) -> Int:
    """Validate records and store six offsets/lengths per record.

    Status is a non-negative record count, or -(one-based FASTQ line number).
    A final newline is optional, matching dnaio's parser.
    """
    var pos = 0
    var records = 0
    while pos < n:
        var name_line = ascii_line_end(data, n, pos)
        if name_line == -2:
            return -(records * 4 + 1)
        if name_line < 0:
            name_line = n
        if pos >= n or data[pos] != UInt8(64):
            return -(records * 4 + 1)
        var seq_start = name_line + 1
        if name_line == n:
            return -(records * 4 + 2)
        var seq_line = ascii_line_end(data, n, seq_start)
        if seq_line < 0:
            return -(records * 4 + 2)
        var plus_start = seq_line + 1
        var plus_line = ascii_line_end(data, n, plus_start)
        if plus_line < 0 or plus_start >= n or data[plus_start] != UInt8(43):
            return -(records * 4 + 3)
        var qual_start = plus_line + 1
        var qual_line = ascii_line_end(data, n, qual_start)
        if qual_line == -2:
            return -(records * 4 + 4)
        if qual_line < 0:
            qual_line = n
        var name_start = pos + 1
        var name_len = trimmed_length(data, name_start, name_line)
        var seq_len = trimmed_length(data, seq_start, seq_line)
        var second_start = plus_start + 1
        var second_len = trimmed_length(data, second_start, plus_line)
        var qual_len = trimmed_length(data, qual_start, qual_line)
        if second_len != 0 and (second_len != name_len or not equal_bytes(data, name_start, second_start, name_len)):
            return -(records * 4 + 3)
        if seq_len != qual_len:
            return -(records * 4 + 4)
        if records < capacity:
            var base = records * 6
            fields[base] = name_start
            fields[base + 1] = name_len
            fields[base + 2] = seq_start
            fields[base + 3] = seq_len
            fields[base + 4] = qual_start
            fields[base + 5] = qual_len
        records += 1
        pos = qual_line + 1
    return records


def copy_bytes(source: BPtr, source_start: Int, width: Int, dest: BPtr, dest_start: Int) -> Int:
    comptime W = simd_width_of[DType.uint8]()
    var i = 0
    while i + W <= width:
        dest.store(dest_start + i, source.load[width=W](source_start + i))
        i += W
    while i < width:
        dest[dest_start + i] = source[source_start + i]
        i += 1
    return dest_start + width


@export("md_scan_fastq")
def md_scan_fastq(data_addr: Int, n: Int, fields_addr: Int, capacity: Int) abi("C") -> Int:
    # The Python bridge always supplies non-null, contiguous buffers. Validate the
    # C ABI inputs before turning integer addresses into non-null Mojo pointers.
    if data_addr == 0 or fields_addr == 0 or n < 0 or capacity < 0:
        return -1
    return scan_fastq(BPtr(unsafe_from_address=data_addr), n,
                      IPtr(unsafe_from_address=fields_addr), capacity)


@export("md_format_fastq")
def md_format_fastq(name_addr: Int, name_len: Int, sequence_addr: Int, sequence_len: Int,
                    qualities_addr: Int, quality_len: Int, two_headers: Int, dest_addr: Int,
                    dest_capacity: Int) abi("C") -> Int:
    if name_addr == 0 or sequence_addr == 0 or qualities_addr == 0 or dest_addr == 0:
        return -2
    if name_len < 0 or sequence_len < 0 or quality_len < 0 or dest_capacity < 0:
        return -2
    if sequence_len != quality_len:
        return -1
    # Check capacity by subtraction so hostile length values cannot overflow a
    # summed required-size calculation before any pointer is dereferenced.
    if dest_capacity < 6:
        return -2
    var remaining = dest_capacity - 6
    if name_len > remaining:
        return -2
    remaining -= name_len
    if sequence_len > remaining:
        return -2
    remaining -= sequence_len
    if quality_len > remaining:
        return -2
    remaining -= quality_len
    if two_headers != 0 and name_len > remaining:
        return -2
    var name = BPtr(unsafe_from_address=name_addr)
    var sequence = BPtr(unsafe_from_address=sequence_addr)
    var qualities = BPtr(unsafe_from_address=qualities_addr)
    var dest = BPtr(unsafe_from_address=dest_addr)
    var p = 0
    dest[p] = UInt8(64)
    p += 1
    p = copy_bytes(name, 0, name_len, dest, p)
    dest[p] = UInt8(10)
    p += 1
    p = copy_bytes(sequence, 0, sequence_len, dest, p)
    dest[p] = UInt8(10)
    p += 1
    dest[p] = UInt8(43)
    p += 1
    if two_headers != 0:
        p = copy_bytes(name, 0, name_len, dest, p)
    dest[p] = UInt8(10)
    p += 1
    p = copy_bytes(qualities, 0, quality_len, dest, p)
    dest[p] = UInt8(10)
    return p + 1
