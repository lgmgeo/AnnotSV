from __future__ import annotations

from typing import TYPE_CHECKING
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context

## - Check if the following chromFa files has been downloaded:
#    - *chromFa.tar.gz
#
## - Check and create if necessary the following file:
#    - 'date'_genomeBuild_chromFa.fasta
#    - 'date'_genomeBuild_chromFa.chromSizes


class FastaValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="GCcontent",
            downloaded=ResolvedFiles(app.config.gcc_dir, "*chromFa.tar.gz"),
            formatted=ResolvedFiles(app.config.gcc_dir, f"{app.config.genome_build}_chromFa.fasta"),
        )

    def update(self):
        raise NotImplementedError()


def gc_content_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
