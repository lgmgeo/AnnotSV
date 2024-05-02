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
import shutil
import sys

from os.path import join as osj

import variantconvert
from test_main import identical_except_date_and_genome, remove_if_exists

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from variantconvert.__main__ import main_convert
from variantconvert.config import change_config


def annotsv_test_routine(
    tmp_path: str, input_file: str, control: str, config: str, config_changes={}
) -> None:
    config_test = osj(tmp_path, os.path.basename(config))
    shutil.copy2(config, config_test)
    if config_changes != {}:
        change_config(config_test, config_changes)

    output_base = os.path.basename(input_file)
    output_name = output_base[: output_base.rfind(".")] + ".vcf"
    annotsv_tester = type(
        "obj",
        (object,),
        {
            "inputFile": input_file,
            "outputFile": osj(tmp_path, output_name),
            "configFile": config_test,
            "verbosity": "debug",
        },
    )
    remove_if_exists(annotsv_tester.outputFile)
    main_convert(annotsv_tester)
    assert os.path.exists(annotsv_tester.outputFile)

    identical_except_date_and_genome([control, annotsv_tester.outputFile])


def test_annotsv_to_vcf_full(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "DECON.results_all.AnnotSV.tsv")
    control = osj(os.path.dirname(__file__), "controls", "annotsv_full", "decon_annotsv_test.vcf")
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_vcf.json")
    config_change = {"GENERAL": {"mode": "full"}}
    annotsv_test_routine(tmp_path, input_file, control, config, config_change)


def test_annotsv_to_vcf_fullandsplit(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "DECON.results_all.AnnotSV.tsv")
    control = osj(
        os.path.dirname(__file__), "controls", "annotsv_fullandsplit", "decon_annotsv_test.vcf"
    )
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_vcf.json")
    config_change = {"GENERAL": {"mode": "full&split"}}
    annotsv_test_routine(tmp_path, input_file, control, config, config_change)


def test_annotsv_to_vcf_combined(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "DECON.results_all.AnnotSV.tsv")
    control = osj(
        os.path.dirname(__file__), "controls", "annotsv_combined", "decon_annotsv_test.vcf"
    )
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_vcf.json")
    annotsv_test_routine(tmp_path, input_file, control, config)


def test_bed_based_annotsv3_to_vcf_full(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "annotsv_from_bed.tsv")
    control = osj(os.path.dirname(__file__), "controls", "annotsv_full", "annotsv3_from_bed.vcf")
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_bed.json")
    config_change = {"GENERAL": {"mode": "full"}}
    annotsv_test_routine(tmp_path, input_file, control, config, config_change)


def test_bed_based_annotsv3_to_vcf_fullandsplit(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "annotsv_from_bed.tsv")
    control = osj(
        os.path.dirname(__file__), "controls", "annotsv_fullandsplit", "annotsv3_from_bed.vcf"
    )
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_bed.json")
    config_change = {"GENERAL": {"mode": "full&split"}}
    annotsv_test_routine(
        tmp_path,
        input_file,
        control,
        config,
        config_change,
    )


def test_bed_based_annotsv3_to_vcf_combined(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "annotsv_from_bed.tsv")
    control = osj(
        os.path.dirname(__file__), "controls", "annotsv_combined", "annotsv3_from_bed.vcf"
    )
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_bed.json")
    annotsv_test_routine(tmp_path, input_file, control, config)


def test_multisample_bed_based_annotsv3_to_vcf_full(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "multisample_from_bed.annotated.tsv")
    control = osj(
        os.path.dirname(__file__), "controls", "annotsv_full", "multisample_annotsv3_from_bed.vcf"
    )
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_bed.json")
    config_change = {"GENERAL": {"mode": "full"}}
    annotsv_test_routine(tmp_path, input_file, control, config, config_change)


def test_multisample_bed_based_annotsv3_to_vcf_fullandsplit(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "multisample_from_bed.annotated.tsv")
    control = osj(
        os.path.dirname(__file__),
        "controls",
        "annotsv_fullandsplit",
        "multisample_annotsv3_from_bed.vcf",
    )
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_bed.json")
    config_change = {"GENERAL": {"mode": "full&split"}}
    annotsv_test_routine(tmp_path, input_file, control, config, config_change)


def test_multisample_bed_based_annotsv3_to_vcf_combined(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "multisample_from_bed.annotated.tsv")
    control = osj(
        os.path.dirname(__file__),
        "controls",
        "annotsv_combined",
        "multisample_annotsv3_from_bed.vcf",
    )
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_bed.json")
    annotsv_test_routine(tmp_path, input_file, control, config)


def test_annotsv_with_wild_types_to_vcf_full(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "annotsv_wt_samples.tsv")
    control = osj(os.path.dirname(__file__), "controls", "annotsv_full", "annotsv_wt_samples.vcf")
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_vcf.json")
    config_change = {"GENERAL": {"mode": "full"}}
    annotsv_test_routine(tmp_path, input_file, control, config, config_change)


def test_annotsv_with_wild_types_to_vcf_fullandsplit(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "annotsv_wt_samples.tsv")
    control = osj(
        os.path.dirname(__file__), "controls", "annotsv_fullandsplit", "annotsv_wt_samples.vcf"
    )
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_vcf.json")
    config_change = {"GENERAL": {"mode": "full&split"}}
    annotsv_test_routine(tmp_path, input_file, control, config, config_change)


def test_annotsv_with_wild_types_to_vcf_full(tmp_path):
    input_file = osj(os.path.dirname(__file__), "data", "annotsv_wt_samples.tsv")
    control = osj(
        os.path.dirname(__file__), "controls", "annotsv_combined", "annotsv_wt_samples.vcf"
    )
    config = osj(variantconvert.__default_config__, "hg19", "annotsv3_from_vcf.json")
    annotsv_test_routine(tmp_path, input_file, control, config)
