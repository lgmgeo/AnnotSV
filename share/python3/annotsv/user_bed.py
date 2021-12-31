from __future__ import annotations

from pathlib import Path
from annotsv.context import Context
from annotsv.constants import install_dir


def check_users_bed_files(app: Context):
    bed_files = list(app.config.user_bed_dir.glob("*/*.bed"))

    if not bed_files:
        app.log.debug(
            f"No user bed annotation in {app.config.user_bed_dir.relative_to(install_dir)}"
        )
    else:
        for bed_file in bed_files:
            if _is_unformatted(bed_file):
                raise NotImplementedError()
            elif bed_file.name.endswith(".formatted.bed"):
                app.log.warning(f"Ignoring intermediary user BED file: {bed_file}")
                continue

            header_file = bed_file.with_suffix(".header.tsv")
            if not header_file.exists():
                raise FileNotFoundError(
                    f"Missing required header file for {bed_file.relative_to(app.config.annotation_root)}: {header_file.name}"
                )

            # Check that the bed and header files have the same number of columns (and more than 3)
            with header_file.open("rt") as fh:
                num_cols_header = len(fh.readline().split("\t"))
            with bed_file.open("rt") as fh:
                num_cols_bed = len(fh.readline().split("\t"))

            if num_cols_header != num_cols_bed:
                app.log.warning(
                    f"{header_file.name} has {num_cols_header} fields, but {num_cols_bed} were "
                    "expected. Header bed file not used."
                )
                header_file.rename(f"{header_file.name}.bad")
                _write_empty_header(bed_file, num_cols_bed)
            elif num_cols_header <= 3:
                app.log.warning(
                    f"{header_file.name} has {num_cols_header} fields. More than 3 are expected. "
                    "Header bed file not used."
                )
                _write_empty_header(bed_file, num_cols_bed)


def _write_empty_header(bed_file: Path, num_cols: int):
    new_header = bed_file.with_suffix(".header.tsv")
    new_header.write_text("\t".join(["" for _ in range(num_cols)]) + "\n")


def _is_unformatted(file: Path):
    return not file.name.endswith(".sorted.bed") and not file.name.endswith(".formatted.bed")


def user_bed_annotation(
    app: Context,
    user_bedfile: Path,
    sv_chrom: str,
    sv_start: int,
    sv_end: int,
):
    ...
