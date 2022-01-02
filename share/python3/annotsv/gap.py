from __future__ import annotations
from typing import TYPE_CHECKING
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class GapValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="gap",
            downloaded=ResolvedFiles(app.config.gap_dir, "Gap.bed"),
            formatted=ResolvedFiles(app.config.gap_dir, "*_Gap.sorted.bed"),
        )

    def update(self):
        raise NotImplementedError()


def gap_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
