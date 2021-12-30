from __future__ import annotations

from annotsv.context import Context


def check_blacklist_file(app: Context):
    downloaded_file = app.config.blacklist_dir / "ENCODEblacklist.bed"
    formatted_files = list(app.config.blacklist_dir.glob("*_ENCODEblacklist.sorted.bed"))
    label = "ENCODE blacklist"

    if app.blacklist_ann is False:
        app.log.debug(f"No {label} data in output columns, ignoring")
    elif not downloaded_file.exists() and not formatted_files:
        app.log.debug(f"No {label} annotation")
        app.blacklist_ann = False
    elif len(formatted_files) > 1:
        app.keep_last_file(label, formatted_files)
    elif not formatted_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new {label} annotation to format")


def blacklist_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
