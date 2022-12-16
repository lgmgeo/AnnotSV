# -*- coding: utf-8 -*-
"""
#for testing:
docker run --rm -ti --entrypoint=bash -v /home1:/home1 tsvconvert
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/variantconvert_project/variantconvert/variantconvert/__main__.py varankBatch -i /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/variantconvert_project/examples/TSV/ -o /home1/BAS/nicaises/Tests/variantconvert_batch/new_BBS_from_varank.vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/variantconvert_project/variantconvert/configs/config_varank.json
"""

from __future__ import division
from __future__ import print_function

import glob
import logging as log
import multiprocessing
import os
import subprocess
import time

from os.path import join as osj

from commons import set_log_level
from converter_factory import ConverterFactory


def conversion_worker(args):
    """
    args are contained in a tuple for ease of use with multiprocessing
    """
    varank_tsv, bcftools, bgzip, tabix, json_config, tmp_dir = args
    log.debug("###varank_tsv: " + varank_tsv)
    factory = ConverterFactory()
    converter = factory.get_converter("varank", "vcf", json_config)

    sample_name = converter.get_sample_name(varank_tsv)
    log.debug("###sample_name: " + sample_name)
    sample_output = osj(tmp_dir, sample_name + "_from_varank.vcf")
    coords_file = osj(os.path.dirname(varank_tsv), "VCF_Coordinates_Conversion.tsv")

    converter.set_coord_conversion_file(coords_file)
    converter.convert(varank_tsv, sample_output)
    del converter  # otherwise they accumulate in memory until the end of pool.map()

    subprocess.run(
        bcftools + " sort " + sample_output + " > " + sample_output + ".sorted.vcf",
        shell=True,
    )
    subprocess.run(bgzip + " " + sample_output + ".sorted.vcf", shell=True)
    subprocess.run(tabix + " -p vcf " + sample_output + ".sorted.vcf.gz", shell=True)


def main_varank_batch(args):
    set_log_level(args.verbosity)
    tmp_dir = osj(os.path.dirname(args.outputFile), ".varanktovcf." + str(time.time()))
    os.makedirs(tmp_dir)
    files_to_convert = glob.glob(
        osj(args.inputVarankDir, "*_allVariants.rankingByVar.tsv")
    )
    if len(files_to_convert) == 0:
        raise ValueError(
            "Expected to find files with pattern '*_allVariants.rankingByVar.tsv' in directory:"
            + args.inputVarankDir
        )

    # Without multiprocessing for easier debugging
    # for varank_tsv in files_to_convert:
    #     myargs = (
    #         varank_tsv,
    #         args.bcftools,
    #         args.bgzip,
    #         args.tabix,
    #         args.configFile,
    #         tmp_dir,
    #     )
    #     conversion_worker(myargs)

    with multiprocessing.Pool(args.ncores) as pool:
        pool.map(
            conversion_worker,
            [
                (
                    varank_tsv,
                    args.bcftools,
                    args.bgzip,
                    args.tabix,
                    args.configFile,
                    tmp_dir,
                )
                for varank_tsv in files_to_convert
            ],
        )

    cmd = (
        args.bcftools
        + " merge -m none "
        + osj(tmp_dir, "*.vcf.gz")
        + " -o "
        + args.outputFile
        + " --threads "
        + str(args.ncores)
    )
    print(cmd)
    subprocess.run(cmd, shell=True)
    # shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    pass
