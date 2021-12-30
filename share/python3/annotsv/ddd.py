from __future__ import annotations

from annotsv.context import Context


def check_ddd_file(app: Context):
    file_downloaded = app.config.extann_dir / "DDD/DDG2P.csv.gz"
    formatted_files = list(app.config.extann_dir.glob("DDD/*_DDG2P.tsv.gz"))

    if not file_downloaded.exists() and not formatted_files:
        app.log.debug("No DDD gene annotation")
    elif len(formatted_files) > 1:
        app.keep_last_file("DDG2P", formatted_files)
    elif not formatted_files:
        update_ddd_file(app)
    else:
        app.log.debug("No new DDD gene annotation to format")


def update_ddd_file(app: Context):
    raise NotImplementedError()


### Haploinsufficiency
# comes from haploinsufficiency.tcl, but moved to here since it's stored under DDD


def check_hi_file(app: Context):
    file_downloaded = app.config.extann_dir / "DDD/HI_Predictions"
    formatted_files = list(app.config.extann_dir.glob("DDD/*_HI.tsv.gz"))

    if not file_downloaded.exists() and not formatted_files:
        app.log.debug("No haploinsufficiency annotation")
    elif len(formatted_files) > 1:
        app.keep_last_file("Haploinsufficiency", formatted_files)
    elif not formatted_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new haploinsufficiency annotation to format")
