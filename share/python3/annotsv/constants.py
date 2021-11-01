import os
import platform
from pathlib import Path

install_dir = Path(__file__).absolute().parents[3]
annotation_dir = install_dir / "share" / "AnnotSV"

version_file = Path(__file__).parents[3] / "VERSION"
annotsv_version = version_file.read_text().strip()

license_blurb = f"""AnnotSV {annotsv_version}

Copyright (C) 2017-2021 GEOFFROY Veronique

Please feel free to contact me for any suggestions or bug reports
email: veronique.geoffroy@inserm.fr

Python version: {platform.python_version()}

Application name used (defined with the "ANNOTSV" environment variable):
{os.getenv("ANNOTSV")}

COMMAND LINE USAGE

       $ANNOTSV/bin/AnnotSV -SVinputFile 'Path of your VCF or BED input file with SV coordinates' >& AnnotSV.log &
"""
