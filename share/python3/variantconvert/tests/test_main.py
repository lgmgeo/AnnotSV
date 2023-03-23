# -*- coding: utf-8 -*-
"""
Test files are stored outside of the git project for confidentiality.
File structure for tests to work properly:

project_folder
    examples
        file_for_tests1.tsv
        file_for_tests2.tsv
    variantconvert
        .git
        tests
        config
        variantconvert
            __main__.py


For testing purposes, GENOME["path"] was changed in configs. Normally it is:
"path": "/home1/data/STARK/databases/genomes/current/hg19.fa"
"""

import logging as log
import os
import pytest
import sys
import typing

from os.path import join as osj

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from variantconvert.__main__ import main_convert


def identical_except_date_and_genome(vcf_list: typing.List[str]):
    """
    Test if all vcf files in the input list are identical, except for 2 fields:
        - VCF creation date, as this will change all the time
        - Genome path, as this part of config can need a local change
    """
    if len(vcf_list) < 2:
        raise RuntimeError("Expected a list containing at least 2 VCF paths")

    VARIABLE_FIELDS = {"##fileDate=", "##reference=", "##InputFile="}

    vcfs_lines = []
    k = 0
    for vcf in vcf_list:
        vcfs_lines.append([])
        with open(vcf, "r") as f:
            for l in f:
                if all([not l.startswith(field) for field in VARIABLE_FIELDS]):
                    vcfs_lines[k].append(l)
        k += 1

    for i in range(len(vcf_list) - 1):
        assert vcfs_lines[i] == vcfs_lines[i + 1]


def remove_if_exists(filepath: str):
    if os.path.exists(filepath):
        os.remove(filepath)


def test_varank_to_vcf(tmp_path):
    varank_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "data",
                "fam01_SAMPLE_VARANK_hg19_allVariants.rankingByGene.tsv",
            ),
            "outputFile": osj(tmp_path, "varank_test.vcf"),
            "inputFormat": "varank",
            "outputFormat": "vcf",
            "configFile": osj(os.path.dirname(__file__), "..", "configs", "hg19", "varank.json"),
            "verbosity": "debug",
            "coordConversionFile": osj(
                os.path.dirname(__file__),
                "data",
                "VCF_Coordinates_Conversion.tsv",
            ),
        },
    )
    remove_if_exists(varank_tester.outputFile)
    main_convert(varank_tester)
    assert os.path.exists(varank_tester.outputFile)

    control = osj(os.path.dirname(__file__), "controls", "varank_test.vcf")
    identical_except_date_and_genome([control, varank_tester.outputFile])


def test_decon_to_vcf(tmp_path):
    decon_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "data",
                "DECON.Design_results_all.txt",
            ),
            "outputFile": osj(tmp_path, "decon_test.vcf"),
            "inputFormat": "tsv",
            "outputFormat": "vcf",
            "configFile": osj(os.path.dirname(__file__), "..", "configs", "hg19", "decon.json"),
            "verbosity": "debug",
        },
    )
    remove_if_exists(decon_tester.outputFile)
    main_convert(decon_tester)
    assert os.path.exists(decon_tester.outputFile)

    control = osj(os.path.dirname(__file__), "controls", "decon_test.vcf")
    identical_except_date_and_genome([control, decon_tester.outputFile])


def test_annotsv_to_vcf(tmp_path):
    annotsv_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "data",
                "DECON.results_all.AnnotSV.tsv",
            ),
            "outputFile": osj(tmp_path, "decon_annotsv_test.vcf"),
            "inputFormat": "annotsv",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__), "..", "configs", "hg19", "annotsv3_from_vcf.json"
            ),
            "verbosity": "debug",
        },
    )
    remove_if_exists(annotsv_tester.outputFile)
    main_convert(annotsv_tester)
    assert os.path.exists(annotsv_tester.outputFile)

    control = osj(os.path.dirname(__file__), "controls", "decon_annotsv_test.vcf")
    identical_except_date_and_genome([control, annotsv_tester.outputFile])


def test_bed_to_vcf(tmp_path):
    bed_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "data",
                "canoes.bed",
            ),
            "outputFile": osj(tmp_path, "canoes_bed.vcf"),
            "inputFormat": "tsv",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__), "..", "configs", "hg19", "canoes_bed.json"
            ),
            "verbosity": "debug",
        },
    )
    remove_if_exists(bed_tester.outputFile)
    main_convert(bed_tester)
    assert os.path.exists(bed_tester.outputFile)

    control = osj(os.path.dirname(__file__), "controls", "canoes_bed.vcf")
    identical_except_date_and_genome([control, bed_tester.outputFile])


def test_breakpoints_to_vcf(tmp_path):
    breakpoints_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "data",
                "star-fusion.fusion_predictions.tsv",
            ),
            "outputFile": osj(tmp_path, "star-fusion.vcf"),
            "inputFormat": "breakpoints",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__), "..", "configs", "hg19", "starfusion.json"
            ),
            "verbosity": "debug",
        },
    )
    remove_if_exists(breakpoints_tester.outputFile)
    main_convert(breakpoints_tester)
    assert os.path.exists(breakpoints_tester.outputFile)

    control = osj(os.path.dirname(__file__), "controls", "star-fusion.vcf")
    identical_except_date_and_genome([control, breakpoints_tester.outputFile])


def test_arriba_breakpoints_to_vcf(tmp_path):
    breakpoints_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "data",
                "arriba.fusions.tsv",
            ),
            "outputFile": osj(tmp_path, "arriba.vcf"),
            "inputFormat": "breakpoints",
            "outputFormat": "vcf",
            "configFile": osj(os.path.dirname(__file__), "..", "configs", "hg19", "arriba.json"),
            "verbosity": "debug",
        },
    )
    remove_if_exists(breakpoints_tester.outputFile)
    main_convert(breakpoints_tester)
    assert os.path.exists(breakpoints_tester.outputFile)

    control = osj(os.path.dirname(__file__), "controls", "arriba.vcf")
    identical_except_date_and_genome([control, breakpoints_tester.outputFile])


def test_bed_based_annotsv3_to_vcf(tmp_path):
    breakpoints_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "data",
                "annotsv_from_bed.tsv",
            ),
            "outputFile": osj(tmp_path, "annotsv3_from_bed.vcf"),
            "inputFormat": "annotsv",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__), "..", "configs", "hg19", "annotsv3_from_bed.json"
            ),
            "verbosity": "debug",
        },
    )
    remove_if_exists(breakpoints_tester.outputFile)
    main_convert(breakpoints_tester)
    assert os.path.exists(breakpoints_tester.outputFile)

    control = osj(os.path.dirname(__file__), "controls", "annotsv3_from_bed.vcf")
    identical_except_date_and_genome([control, breakpoints_tester.outputFile])


def test_multisample_bed_based_annotsv3_to_vcf(tmp_path):
    breakpoints_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "data",
                "multisample_from_bed.annotated.tsv",
            ),
            "outputFile": osj(tmp_path, "multisample_annotsv3_from_bed.vcf"),
            "inputFormat": "annotsv",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__), "..", "configs", "hg19", "annotsv3_from_bed.json"
            ),
            "verbosity": "debug",
        },
    )
    remove_if_exists(breakpoints_tester.outputFile)
    main_convert(breakpoints_tester)
    assert os.path.exists(breakpoints_tester.outputFile)

    control = osj(os.path.dirname(__file__), "controls", "multisample_annotsv3_from_bed.vcf")
    identical_except_date_and_genome([control, breakpoints_tester.outputFile])


def test_bedpe_to_vcf(tmp_path):
    bed_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "data",
                "chromothripsis.bedpe",
            ),
            "outputFile": osj(tmp_path, "chromo.vcf"),
            "inputFormat": "bedpe",
            "outputFormat": "vcf",
            "configFile": osj(os.path.dirname(__file__), "..", "configs", "hg19", "bedpe.json"),
            "verbosity": "debug",
        },
    )
    remove_if_exists(bed_tester.outputFile)
    main_convert(bed_tester)
    assert os.path.exists(bed_tester.outputFile)

    control = osj(os.path.dirname(__file__), "controls", "chromothripsis.vcf")
    identical_except_date_and_genome([control, bed_tester.outputFile])
