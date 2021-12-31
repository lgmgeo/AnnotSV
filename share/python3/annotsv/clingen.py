from __future__ import annotations

import typing

from annotsv.schemas import AnnotationValidator, ResolvedFiles

if typing.TYPE_CHECKING:
    from annotsv.context import Context


class ClingenValidator(AnnotationValidator):
    def __init__(self, app: Context):
        download_rf = ResolvedFiles(app.config.clingen_dir, "ClinGen_gene_curation_list_*.tsv")
        formatted_rf = ResolvedFiles(app.config.clingen_dir, "*_ClinGenAnnotations.tsv.gz")
        super().__init__(
            app,
            label="ClinGen",
            downloaded=download_rf,
            formatted=formatted_rf,
        )

    def update(self):
        raise NotImplementedError()
