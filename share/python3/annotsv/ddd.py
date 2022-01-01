from __future__ import annotations

from typing import TYPE_CHECKING
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class DDDValidator(AnnotationValidator):
    def __init__(self, app: Context):
        label = "DDG2P"
        downloaded_rf = ResolvedFiles(app.config.ddd_dir, "DDG2P.csv.gz")
        formatted_rf = ResolvedFiles(app.config.ddd_dir, "*_DDG2P.tsv.gz")
        super().__init__(app, label=label, downloaded=downloaded_rf, formatted=formatted_rf)

    def update(self):
        raise NotImplementedError()


### Haploinsufficiency
# comes from haploinsufficiency.tcl, but moved to here since it's stored under DDD
class HaploinsufficiencyValidator(AnnotationValidator):
    def __init__(self, app: Context):
        label = "Haploinsufficiency"
        downloaded_rf = ResolvedFiles(app.config.ddd_dir, "HI_Predictions")
        formatted_rf = ResolvedFiles(app.config.ddd_dir, "*_HI.tsv.gz")
        super().__init__(app, label=label, downloaded=downloaded_rf, formatted=formatted_rf)

    def update(self):
        raise NotImplementedError()
