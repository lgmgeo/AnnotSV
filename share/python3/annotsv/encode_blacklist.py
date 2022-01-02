from __future__ import annotations
from typing import TYPE_CHECKING
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class BlacklistValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="ENCODE blacklist",
            downloaded=ResolvedFiles(app.config.blacklist_dir, "ENCODEblacklist.bed"),
            formatted=ResolvedFiles(app.config.blacklist_dir, "*_ENCODEblacklist.sorted.bed"),
        )

    def update(self):
        raise NotImplementedError()


def blacklist_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
