from __future__ import annotations

from annotsv.context import Context


def check_loeuf_file(app: Context):
    downloaded_files = list(app.config.extann_dir.glob("gnomAD/gnomad.*.lof_metrics.by_gene.txt"))
    formatted_files = list(
        app.config.extann_dir.glob("gnomAD/*_gnomAD.LOEUF.pLI.annotations.tsv.gz")
    )

    if not downloaded_files and not formatted_files:
        app.log.debug("No LOEUF annotation")
    elif len(formatted_files) > 1:
        app.keep_last_file("LOEUF", formatted_files)
    elif not formatted_files:
        raise NotImplementedError()
    else:
        app.log.debug("No new LOEUF annotation to format")
