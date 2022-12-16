# -*- coding: utf-8 -*-
"""
@Goal: Expand Celine Besnard's script with infinite conversion abilities between vcf and various other formats
@Author: Samuel Nicaise
@Date: 23/11/2021

Prerequisites: pandas, pyfaidx (https://github.com/mdshw5/pyfaidx))

Usage examples:
#TSV (Decon) to VCF
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/file_converter.py -i /home1/BAS/nicaises/Tests/deconconverter/200514_NB551027_0724_AHTWHHAFXY.DECON_results_all.txt -o vcf_from_decon.vcf -fi tsv -fo vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/fileconversion/config_decon.json

#AnnotSV3 to VCF
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/file_converter.py -i /home1/BAS/nicaises/Tests/deconconverter/DECoN.20211207-183955_results_both.tsv -o /home1/BAS/nicaises/Tests/deconconverter/decon__annotsv3.vcf -fi annotsv -fo vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/fileconversion/config_annotsv3.json

#Canoes BED to VCF
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/file_converter.py -i canoes_calling.bed -o vcf_from_canoes_bed.vcf -fi tsv -fo vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/fileconversion/config_canoes_bed.json

----------------------------
Configfile guidelines (JSON)
----------------------------
1)[GENERAL] has 3 important fields.
	#source format name: will show up in VCF meta fields
	#skip_rows: how many rows to skip before we reach indexes.
	This script cannot handle a tsv with unnamed columns (beds are fine)
	#unique_variant_id: useful in multisample files. List the
	columns that are needed to uniquely identify a variant.
2) [VCF_COLUMNS] describe the columns that will go in your VCF
	key: column name in VCF ; value: column name in source format
3) [COLUMNS_DESCRIPTION] describe the tsv columns
	Type and Description fields will be used in the VCF header
4) read HelperFunctions docstring

If you need a place to store variables unrelated to the vcf file (e.g number of CPUs) put them in [GENERAL]

#TODO: add argument mode to change config files (particularly genome)
#TODO: add argument mode to deal with an entire folder of varank files (or varank files in general)
#TODO: refactor with a Variant class and a VCF class
"""
from __future__ import division
from __future__ import print_function

import argparse
import os
import sys

from os.path import join as osj

sys.path.append(os.path.join(os.path.dirname(__file__), "."))
from commons import set_log_level
from converter_factory import ConverterFactory
from varank_batch import main_varank_batch


def main_convert(args):
    set_log_level(args.verbosity)
    if args.inputFormat.lower() == "decon":
        raise ValueError(
            "DECON is handled as a TSV conversion. Use 'tsv' as input format"
        )

    factory = ConverterFactory()
    converter = factory.get_converter(
        args.inputFormat.lower(), args.outputFormat.lower(), args.configFile
    )

    if args.inputFormat == "varank":
        if args.coordConversionFile == "":
            raise ValueError(
                "Converting from a Varank file requires setting the --coordConversionFile argument to an existing file"
            )
        if not os.path.exists(args.coordConversionFile):
            raise ValueError(
                "coordConversionFile does not exist:" + args.coordConversionFile
            )
        converter.set_coord_conversion_file(args.coordConversionFile)

    converter.convert(args.inputFile, args.outputFile)


def main():
    parser = argparse.ArgumentParser(prog="variantconvert")
    subparsers = parser.add_subparsers(help="sub-command help")
    parser_convert = subparsers.add_parser(
        "convert", help="convert a file containing genomic variants to an other format"
    )
    parser_convert.add_argument(
        "-i", "--inputFile", type=str, required=True, help="Input file"
    )
    parser_convert.add_argument(
        "-o", "--outputFile", type=str, required=True, help="Output file"
    )
    parser_convert.add_argument(
        "-fi", "--inputFormat", type=str, required=True, help="Input file format"
    )
    parser_convert.add_argument(
        "-fo", "--outputFormat", type=str, required=True, help="Output file format"
    )
    parser_convert.add_argument(
        "-c",
        "--configFile",
        type=str,
        required=True,
        help="JSON config file describing columns. See script's docstring.",
    )
    parser_convert.add_argument(
        "-cc",
        "--coordConversionFile",
        type=str,
        default="",
        help="Varank coordinate conversion file (only useful if inputFormat=varank)",
    )

    parser_batch = subparsers.add_parser(
        "varankBatch", help="convert an entire folder of Varank files"
    )
    parser_batch.add_argument(
        "-i",
        "--inputVarankDir",
        type=str,
        required=True,
        help="Input directory containing Varank TSV files and VCF_Coordinates_Conversion.tsv",
    )
    parser_batch.add_argument(
        "-o", "--outputFile", type=str, required=True, help="Output file"
    )
    parser_batch.add_argument(
        "-c",
        "--configFile",
        type=str,
        required=True,
        help="JSON config file describing columns. See script's docstring.",
    )
    parser_batch.add_argument(
        "-n",
        "--ncores",
        type=int,
        default=6,
        help="Number of cores for multiprocessing",
    )
    parser_batch.add_argument(
        "-bc",
        "--bcftools",
        type=str,
        default="bcftools",
        help="path to bcftools executable [default : 'bcftools']",
    )
    parser_batch.add_argument(
        "-bg",
        "--bgzip",
        type=str,
        default="bgzip",
        help="path to bgzip executable [default: 'bgzip']",
    )
    parser_batch.add_argument(
        "-ta",
        "--tabix",
        type=str,
        default="tabix",
        help="path to tabix executable [default: 'tabix']",
    )

    parser_config = subparsers.add_parser(
        "config", help="change variables in config files [under construction]"
    )
    parser_config.add_argument("-g", "--genome", type=str, help="genome path")
    parser_config.add_argument(
        "-c",
        "--configFile",
        type=str,
        required=True,
        help="JSON config file describing columns. See script's docstring.",
    )

    for myparser in (parser_convert, parser_batch, parser_config):
        myparser.add_argument(
            "-v", "--verbosity", type=str, default="info", help="Verbosity level"
        )

    args = parser.parse_args()
    if "verbosity" not in args:
        parser.print_help()
    elif "inputVarankDir" in args:
        main_varank_batch(args)
    else:
        main_convert(args)


if __name__ == "__main__":
    main()