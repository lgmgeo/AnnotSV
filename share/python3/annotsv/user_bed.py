from __future__ import annotations

from typing import TYPE_CHECKING
from pathlib import Path
from annotsv.constants import install_dir
from annotsv.schemas import AnnotationValidator, ResolvedFiles
from annotsv.bed import check_bed, write_empty_header

if TYPE_CHECKING:
    from annotsv.context import Context


class UserBEDValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="user BED",
            downloaded=ResolvedFiles(app.config.user_bed_dir, "*/*.bed"),
            formatted=ResolvedFiles(app.config.user_bed_dir, "*/*.formatted.sorted.bed"),
        )

    def downloaded_exist(self):
        return super().downloaded_exist() and self.get_downloaded()

    def get_downloaded(self):
        plist = super().get_downloaded()[0]
        if isinstance(plist, Path):
            plist = [plist]
        return [p for p in plist if _is_unformatted(p)]

    def check(self):
        success = super().check()
        if success:
            fpaths = self.get_formatted()[0]
            if isinstance(fpaths, Path):
                fpaths = [fpaths]
            for bed_file in fpaths:
                header_file = bed_file.with_suffix(".header.tsv")
                if not header_file.exists():
                    check_bed(self._app, bed_file)

                # Check that the bed and header files have the same number of columns (and more than 3)
                with header_file.open("rt") as fh:
                    num_cols_header = len(fh.readline().split("\t"))
                with bed_file.open("rt") as fh:
                    num_cols_bed = len(fh.readline().split("\t"))

                if num_cols_header != num_cols_bed:
                    self._app.log.warning(
                        f"{header_file.name} has {num_cols_header} fields, but {num_cols_bed} were "
                        "expected. Header bed file not used."
                    )
                    header_file.rename(f"{header_file.name}.bad")
                    write_empty_header(bed_file, num_cols_bed)
                elif num_cols_header <= 3:
                    self._app.log.warning(
                        f"{header_file.name} has {num_cols_header} fields. More than 3 are expected. "
                        "Header bed file not used."
                    )
                    write_empty_header(bed_file, num_cols_bed)
        else:
            self._app.log.debug(
                f"No user bed annotation in {self._app.config.user_bed_dir.relative_to(install_dir)}"
            )

    def update(self):
        for bed_file in self.get_downloaded():
            check_bed(self._app, bed_file)


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
                write_empty_header(bed_file, num_cols_bed)
            elif num_cols_header <= 3:
                app.log.warning(
                    f"{header_file.name} has {num_cols_header} fields. More than 3 are expected. "
                    "Header bed file not used."
                )
                write_empty_header(bed_file, num_cols_bed)


def _is_unformatted(file: Path):
    return (
        file.name.endswith(".bed")
        and not file.name.endswith(".sorted.bed")
        and not file.name.endswith(".formatted.bed")
    )


def user_bed_annotation(
    app: Context,
    user_bedfile: Path,
    sv_chrom: str,
    sv_start: int,
    sv_end: int,
):
    ...
