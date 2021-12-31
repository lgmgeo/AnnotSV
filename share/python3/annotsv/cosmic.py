from __future__ import annotations

from annotsv.context import Context


def check_cosmic_file(app: Context):
    downloaded_file = app.config.cosmic_dir / "CosmicCompleteCNA.tsv.gz"
    formatted_file = app.config.cosmic_dir / f"CosmicCompleteCNA_{app.config.genome_build}.bed"
    label = "COSMIC"

    if formatted_file.exists():
        app.log.debug(f"Enabling {label} annotation")
        app.cosmic_ann = True
    elif not downloaded_file.exists():
        app.log.debug(f"No {label} annotation")
    else:
        raise NotImplementedError()
