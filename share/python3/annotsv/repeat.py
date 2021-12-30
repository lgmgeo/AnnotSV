from __future__ import annotations

from annotsv.context import Context


def check_repeat_file(app: Context):
    downloaded_file = app.config.repeat_dir / "Repeat.bed"
    formatted_files = list(app.config.repeat_dir.glob("*_Repeat.sorted.bed"))
    label = "Repeat"

    if not downloaded_file.exists() and not formatted_files:
        app.log.debug(f"No {label} annotation")
        app.repeat_ann = False
    elif app.repeat_ann is False:
        app.log.debug(f"No {label} data in output columns, ignoring")
    elif len(formatted_files) > 1:
        app.keep_last_file(label, formatted_files)
    elif not formatted_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new {label} annotation to format")


def repeat_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
