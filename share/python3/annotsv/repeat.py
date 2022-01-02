from __future__ import annotations

from typing import TYPE_CHECKING
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class RepeatValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="Repeat",
            downloaded=ResolvedFiles(app.config.repeat_dir, "Repeat.bed"),
            formatted=ResolvedFiles(app.config.repeat_dir, "*_Repeat.sorted.bed"),
        )

    def update(self):
        raise NotImplementedError()


def repeat_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
