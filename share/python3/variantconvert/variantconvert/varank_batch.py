# -*- coding: utf-8 -*-
"""
@Author: Samuel Nicaise
@Version: v1.0.0

#for testing:
docker run --rm -ti --entrypoint=bash -v /home1:/home1 tsvconvert
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/variantconvert_project/variantconvert/variantconvert/__main__.py varankBatch -i /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/variantconvert_project/examples/TSV/ -o /home1/BAS/nicaises/Tests/variantconvert_batch/new_BBS_from_varank.vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/variantconvert_project/variantconvert/configs/config_varank.json
"""

import glob
import logging as log
import multiprocessing
import os
import tqdm
import time

from os.path import join as osj

from commons import run_shell, set_log_level
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

    run_shell(bcftools + " sort " + sample_output + " > " + sample_output + ".sorted.vcf")
    run_shell(bgzip + " " + sample_output + ".sorted.vcf")
    run_shell(tabix + " -p vcf " + sample_output + ".sorted.vcf.gz")


def main_varank_batch(args):
    set_log_level(args.verbosity)
    tmp_dir = osj(os.path.dirname(args.outputFile), ".varanktovcf." + str(time.time()))
    os.makedirs(tmp_dir)
    files_to_convert = glob.glob(osj(args.inputVarankDir, "*_allVariants.rankingByVar.tsv"))
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

    worker_args = [
        (
            varank_tsv,
            args.bcftools,
            args.bgzip,
            args.tabix,
            args.configFile,
            tmp_dir,
        )
        for varank_tsv in files_to_convert
    ]

    with multiprocessing.Pool(args.ncores) as pool:
        if log.root.level <= 20:
            with tqdm.tqdm(total=len(worker_args), desc="Converting files") as pbar:
                for _ in pool.imap_unordered(conversion_worker, worker_args):
                    pbar.update()
        else:
            # don't show progress bar if log_level >= WARNING
            pool.map(conversion_worker, worker_args)

    if len(files_to_convert) > args.max_merged:
        # Issue: bcftools is not able to merge too many files at once (depends on local system)
        # Solution: merge vcf per subsets of args.max_merged files at once
        # Then merge the merged vcf subsets
        commands = [
            f"ls {osj(tmp_dir, '*.vcf.gz')} | split -l {args.max_merged} - {osj(tmp_dir, 'subset_vcfs')}",
            f"for s in {osj(tmp_dir, 'subset_vcfs*')}; do {args.bcftools} merge -m none --threads {args.ncores} -l $s -o {osj(tmp_dir, 'merge.$(basename $s).vcf')}; {args.bgzip} {osj(tmp_dir, 'merge.$(basename $s).vcf')}; {args.tabix} -p vcf {osj(tmp_dir, 'merge.$(basename $s).vcf.gz')}; done",
            f"ls {osj(tmp_dir, 'merge.*.vcf.gz')} > {osj(tmp_dir, 'merge_list.txt')}",
            f"{args.bcftools} merge --threads {args.ncores} -l {osj(tmp_dir, 'merge_list.txt')} -m none -o {args.outputFile}",
        ]
        for cmd in commands:
            log.info("Merging converted files...")
            log.debug(cmd)
            run_shell(cmd)

    else:
        # merge everything at once
        cmd = f"{args.bcftools} merge -m none {osj(tmp_dir, '*.vcf.gz')} -o {args.outputFile} --threads {args.ncores}"
        log.info("Merging converted files...")
        log.debug(cmd)
        run_shell(cmd)
    # shutil.rmtree(tmp_dir)
