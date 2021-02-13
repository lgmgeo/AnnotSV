from pathlib import Path

install_dir = Path(__file__).absolute().parents[3]
annotation_dir = install_dir / "share" / "AnnotSV"

license_blurb = """AnnotSV {3.0.5}

Copyright (C) 2017-2021 GEOFFROY Veronique

Please feel free to contact me for any suggestions or bug reports
email: veronique.geoffroy@inserm.fr

Tcl/Tk version: {8.6}

Application name used (defined with the "ANNOTSV" environment variable):
{/home/ruphos/playground/AnnotSV}



COMMAND LINE USAGE

       $ANNOTSV/bin/AnnotSV -SVinputFile 'Path of your VCF or BED input file with SV coordinates' >& AnnotSV.log &
"""
