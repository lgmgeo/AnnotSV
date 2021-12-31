from __future__ import annotations
from typing import List

from annotsv.context import Context


def check_tad_files(app: Context):
    downloaded_files = list(app.config.tad_dir.glob("ENC*.bed"))
    formatted_files = list(app.config.tad_dir.glob("*_boundariesTAD.sorted.bed"))
    label = "TAD"

    if app.tad_ann is False:
        app.log.debug(f"No {label} data in output columns, ignoring")
    elif not downloaded_files and not formatted_files:
        app.log.debug(f"No {label} annotation")
        app.tad_ann = False
    elif len(formatted_files) > 1:
        app.keep_last_file(label, formatted_files)
    elif not formatted_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new {label} annotation to format")


def tad_annotation(app: Context, sv_chrom: str, sv_start: int, sv_end: int, ilist: List):
    ...
