from __future__ import annotations
from typing import TYPE_CHECKING
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class LoeufValidator(AnnotationValidator):
    def __init__(self, app: Context):
        # , *, label: str, downloaded: ResolvedFiles, formatted: ResolvedFiles, extra_downloaded: List[ResolvedFiles] = None, extra_formatted: List[ResolvedFiles] = None):
        label = "LOEUF"
        downloaded_rf = ResolvedFiles(
            app.config.gnomad_dir, "gnomAD/gnomad.*.lof_metrics.by_gene.txt"
        )
        formatted_rf = ResolvedFiles(app.config.gnomad_dir, "*_gnomAD.LOEUF.pLI.annotations.tsv.gz")
        super().__init__(app, label=label, downloaded=downloaded_rf, formatted=formatted_rf)

    def update(self) -> bool:
        raise NotImplementedError()
