class UnknownFileFormat(Exception):
    """The file format could not be automatically detected."""


class FileFormatError(Exception):
    format = "sequence"

    def __init__(self, msg: str, line: int | None):
        super().__init__(msg, line)
        self.message, self.line = msg, line

    def __str__(self) -> str:
        where = "unknown line" if self.line is None else f"line {self.line + 1}"
        return f"Error in {self.format} file at {where}: {self.message}"


class FastqFormatError(FileFormatError):
    format = "FASTQ"
