from collections import defaultdict
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Optional, Sequence, TYPE_CHECKING, Set

from annotsv.constants import VALID_CHROM

# if TYPE_CHECKING:
from annotsv.context import Context

CHROM = 0
START = 1
END = 2
ANNOTATION = 3


class BedReader:
    def __init__(self, file: Path):
        self.file = file
        if self.file.suffix == ".gz":
            open_func = gzip.open
        else:
            open_func = open
        with open_func(self.file, "rt") as fh:
            for line in fh:
                if not line.strip():
                    continue
                if line.startswith("#"):
                    self.header = [c.strip() for c in line.split("\t")]
                else:
                    self.header = []
                break
        self.fh = open_func(self.file, "rt")
        if self.header:
            next(self.fh)

    def __iter__(self):
        return self

    def __next__(self):
        line = ""
        while not line.strip():
            line = next(self.fh)
        cols = line.rstrip("\n\r").split("\t")
        if len(cols) > 3:
            anno = cols[3:]
        else:
            anno = None
        return BedRow(cols[CHROM], int(cols[START]), int(cols[END]), anno)


@dataclass(frozen=True)
class BedRow:
    __slots__ = ["chrom", "start", "end", "annotation"]
    chrom: str
    start: int
    end: int
    annotation: Optional[Sequence[str]]

    def __str__(self) -> str:
        line_str = self.short_str
        if self.annotation:
            line_str = "\t".join([line_str, *self.annotation])
        return line_str

    @property
    def length(self):
        row_len = 3
        if self.annotation:
            row_len += len(self.annotation)
        return row_len

    @property
    def short_str(self):
        return f"{self.chrom}\t{self.start}\t{self.end}"

    def replace(
        self,
        *,
        chrom: str = None,
        start: int = None,
        end: int = None,
        annotation: Sequence[str] = None,
        empty_annotation: bool = False,
    ):
        new_chrom = self.chrom if chrom is None else chrom
        new_start = self.start if start is None else start
        new_end = self.end if end is None else end
        new_anno = self.annotation if annotation is None else annotation
        if empty_annotation:
            new_anno = None
        return self.__class__(new_chrom, new_start, new_end, new_anno)


def check_bed(app: Context, bedfile: Path, new_dir: Path = None, merge_overlap: bool = False):
    if new_dir is None:
        new_dir = bedfile.parent
    fbed = new_dir / bedfile.with_suffix(".formatted.bed").name
    header_file = new_dir / bedfile.with_suffix(".header.tsv").name

    if fbed.exists():
        return fbed

    by_chrom: DefaultDict[str, Set[BedRow]] = defaultdict(set)
    rdr = BedReader(bedfile)
    header = rdr.header

    bad_format = False
    row_length = -1
    for row in rdr:
        if row.chrom.startswith("chr"):
            row = row.replace(chrom=row.chrom[3:])

        if row_length < 0:
            row_length = row.length
        elif row.length != row_length:
            app.log.warning(
                f"{bedfile}: Unexpected file format, not the same length for each line. Not using the associated annotations."
            )
            bad_format = True
            break

        by_chrom[row.chrom].add(row)

    if bad_format:
        use_ann = False
    elif row_length > 3:
        use_ann = True
    else:
        use_ann = False

    if not header_file.exists():
        if rdr.header:
            header_file.write_text(f"{header}\n")
        elif use_ann:
            header_file.write_text("\t".join("" for _ in range(row_length)) + "\n")

    with fbed.open("wt") as fh:
        for chrom in VALID_CHROM:
            if chrom not in by_chrom:
                continue
            lines = sorted(by_chrom[chrom], key=lambda x: x.end)
            lines.sort(key=lambda x: x.start)

            if merge_overlap:
                new_lines = [lines[0]]
                for row in lines[1:]:
                    prev_row = lines[-1]
                    if row.start < prev_row.end:
                        if row.end > prev_row.end:
                            new_lines[-1] = prev_row.replace(end=row.end)
                lines = new_lines

            for line in lines:
                if not use_ann:
                    line = line.replace(empty_annotation=True)
                fh.write(f"{line}\n")

    return fbed


def write_empty_header(bed_file: Path, num_cols: int):
    new_header = bed_file.with_suffix(".header.tsv")
    new_header.write_text("\t".join(["" for _ in range(num_cols)]) + "\n")
