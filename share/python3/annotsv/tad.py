from __future__ import annotations

from typing import TYPE_CHECKING, List

from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class TadValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="TAD",
            downloaded=ResolvedFiles(app.config.tad_dir, "ENC*.bed"),
            formatted=ResolvedFiles(app.config.tad_dir, "*_boundariesTAD.sorted.bed"),
        )

    def update(self):
        raise NotImplementedError()


def tad_annotation(app: Context, sv_chrom: str, sv_start: int, sv_end: int, ilist: List):
    ...
