from __future__ import annotations

from annotsv.context import Context

## - Check if the following chromFa files has been downloaded:
#    - *chromFa.tar.gz
#
## - Check and create if necessary the following file:
#    - 'date'_genomeBuild_chromFa.fasta
#    - 'date'_genomeBuild_chromFa.chromSizes
def check_fasta_files(app: Context):
    downloaded_files = list(app.config.gcc_dir.glob("*chromFa.tar.gz"))
    formatted_file = app.config.gcc_dir / f"{app.config.genome_build}_chromFa.fasta"

    if not downloaded_files and not formatted_file.exists():
        app.log.debug("No GCcontent annotation")
        app.gccontent_ann = False
    elif app.gccontent_ann is False:
        app.log.debug("No GCcontent in output columns, ignoring")
    elif not formatted_file.exists():
        update_fasta_files(app)
    else:
        app.log.debug("No new GCcontent annotation to format")


def update_fasta_files(app: Context):
    raise NotImplementedError()


def gc_content_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
