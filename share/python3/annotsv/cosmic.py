from __future__ import annotations

from typing import TYPE_CHECKING

from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class CosmicValidator(AnnotationValidator):
    def __init__(self, app: Context):
        rf_dir = app.config.cosmic_dir
        super().__init__(
            app,
            label="COSMIC",
            downloaded=ResolvedFiles(rf_dir, "CosmicCompleteCNA.tsv.gz"),
            formatted=ResolvedFiles(rf_dir, f"CosmicCompleteCNA_{app.config.genome_build}.bed"),
        )

    def update(self):
        raise NotImplementedError()
