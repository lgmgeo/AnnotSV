# -*- coding: utf-8 -*-
"""
@Goal: Expand Celine Besnard's script with infinite conversion abilities between vcf and various other formats
@Author: Samuel Nicaise
@Date: 23/11/2021
@Version: v2.0.1

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
"""

import argparse
import json
import logging as log
import os
import variantconvert
import sys

from variantconvert.commons import set_log_level
from variantconvert.config import get_full_path, main_config
from variantconvert.converter_factory import ConverterFactory
from variantconvert.initialize import main_init
from variantconvert.varank_batch import main_varank_batch


def main_convert(args):
    set_log_level(args.verbosity)

    config_path = get_full_path(args.configFile)
    log.debug(f"Using config: {config_path}")
    with open(config_path, "r") as f:
        config = json.load(f)
    input_format = config["GENERAL"].get("input_format", "").lower()
    output_format = config["GENERAL"].get("output_format", "").lower()

    if input_format == "decon":
        raise ValueError("DECON is handled as a TSV conversion. Use 'tsv' as input format")

    factory = ConverterFactory()
    converter = factory.get_converter(input_format, output_format, config_path)

    if input_format == "varank":
        if args.coordConversionFile == "":
            raise ValueError(
                "Converting from a Varank file requires setting the --coordConversionFile argument to an existing file"
            )
        if not os.path.exists(args.coordConversionFile):
            raise ValueError("coordConversionFile does not exist:" + args.coordConversionFile)
        converter.set_coord_conversion_file(args.coordConversionFile)

    converter.convert(args.inputFile, args.outputFile)


def main():
    parser = argparse.ArgumentParser(prog="variantconvert")
    parser.add_argument(
        "--version", action="version", version=f"{parser.prog} {variantconvert.__version__}"
    )
    subparsers = parser.add_subparsers(help="sub-command help")

    parser_convert = subparsers.add_parser(
        "convert",
        help="Convert a file containing genomic variants to an other format",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    parser_convert.set_defaults(subparser="convert")
    parser_convert.add_argument("-i", "--inputFile", type=str, required=True, help="Input file")
    parser_convert.add_argument("-o", "--outputFile", type=str, required=True, help="Output file")
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
        "varankBatch",
        help="Convert an entire folder of Varank files",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    parser_batch.set_defaults(subparser="varankBatch")
    parser_batch.add_argument(
        "-i",
        "--inputVarankDir",
        type=str,
        required=True,
        help="Input directory containing Varank TSV files and VCF_Coordinates_Conversion.tsv",
    )
    parser_batch.add_argument("-o", "--outputFile", type=str, required=True, help="Output file")
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
        help="Number of cores for multiprocessing [default: 6]",
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
    parser_batch.add_argument(
        "-mm",
        "--max_merged",
        type=int,
        default=200,
        help="maximum files merged at once by bcftools [default: 200]",
    )

    parser_config = subparsers.add_parser(
        "config",
        help="Change variables in config files",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    parser_config.add_argument(
        "-c",
        "--configFiles",
        type=str,
        default="<default user dir>",
        help="Config file(s) on which changes are applied. Wildcards are accepted. [default: '<script_dir>/configs/*']",
        nargs="*",
    )
    parser_config.set_defaults(subparser="config")
    parser_config.add_argument(
        "-s",
        "--set",
        type=str,
        help="List of variables to set. Example: --set GENOME.assembly=hg19 GENOME.path=/path/hg19.fa",
        nargs="*",
    )
    parser_config.add_argument(
        "--fill_genome_header",
        action="store_true",
        help="Fill the GENOME['vcf_header'] field based on GENOME['path']",
    )
    parser_config.add_argument(
        "--show_config_dir", action="store_true", help="Prints current local config directory"
    )

    parser_init = subparsers.add_parser(
        "init",
        help="Use this first after installing variantconvert. It copies default configs in your user directory so you can edit them",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    parser_init.set_defaults(subparser="init")
    parser_init.add_argument(
        "-d",
        "--dir",
        type=str,
        help="By default, your standard user config directory is detected and INIT will write there. This option allows you to choose another directory. You then have to set yourself an environment variable VARIANTCONVERT_CONFIG containing this path.",
    )

    for myparser in (parser_convert, parser_batch, parser_config, parser_init):
        myparser.add_argument("-v", "--verbosity", type=str, default="info", help="Verbosity level")

    args = parser.parse_args()

    if not hasattr(args, "subparser"):
        parser.print_help()
    else:
        set_log_level(args.verbosity)
        log.debug(f"Args: {str(args)}")
        if args.subparser == "convert":
            log.info(f"running variantconvert {variantconvert.__version__}")
            main_convert(args)
            log.info("variantconvert finished.")

        elif args.subparser == "varankBatch":
            log.info(f"running variantconvert {variantconvert.__version__} (batch mode)")
            main_varank_batch(args)
            log.info("variantconvert finished.")

        elif args.subparser == "config":
            if not (args.set or args.fill_genome_header or args.show_config_dir):
                raise parser_config.error(
                    "one flag is required: --set or --fill_genome_header or --show_config_dir"
                )
            main_config(args)

        elif args.subparser == "init":
            main_init(args)


if __name__ == "__main__":
    main()
