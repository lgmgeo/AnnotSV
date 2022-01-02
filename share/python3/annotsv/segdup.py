from __future__ import annotations

from typing import TYPE_CHECKING

from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class SegDupValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="SegDup",
            downloaded=ResolvedFiles(app.config.segdup_dir, "SegDup.bed"),
            formatted=ResolvedFiles(app.config.segdup_dir, "*_SegDup.sorted.bed"),
        )

    def update(self):
        raise NotImplementedError()


def segdup_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
