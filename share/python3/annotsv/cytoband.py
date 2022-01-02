from __future__ import annotations
from typing import TYPE_CHECKING
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class CytobandValidator(AnnotationValidator):
    def __init__(self, app: Context):
        pattern_stem = f"cytoBand_{app.config.genome_build}"
        super().__init__(
            app,
            label="cytoBand",
            downloaded=ResolvedFiles(app.config.cytoband_dir, f"{pattern_stem}.bed"),
            formatted=ResolvedFiles(
                app.config.cytoband_dir, f"{pattern_stem}.formatted.sorted.bed"
            ),
            extra_formatted=[ResolvedFiles(app.config.cytoband_dir, f"{pattern_stem}.header.tsv")],
        )

    def update(self):
        raise NotImplementedError()


def cytoband_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
