from __future__ import annotations

from annotsv.context import Context


def check_clingen_file(app: Context):
    downloaded_files = list(app.config.clingen_dir.glob("ClinGen_gene_curation_list_*.tsv"))
    formatted_files = list(app.config.clingen_dir.glob("*_ClinGenAnnotations.tsv.gz"))

    if not downloaded_files and not formatted_files:
        app.log.debug("No ClinGen annotation")
    elif len(formatted_files) > 1:
        app.keep_last_file("ClinGen", formatted_files)
    elif not formatted_files:
        raise NotImplementedError()
    else:
        app.log.debug("No new ClinGen annotation to format")
