# -*- coding: utf-8 -*-
"""
@Author: Samuel Nicaise
@Version: v2.0.1
"""

import logging as log
import os
import subprocess
import time

from functools import lru_cache
from pyfaidx import Fasta

import variantconvert


def set_log_level(verbosity):
    verbosity = verbosity.lower()
    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
    }
    if verbosity not in configs.keys():
        raise ValueError(
            f"Unknown verbosity level: {verbosity}\nPlease use any in: {configs.keys()}"
        )
    log.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )


def run_shell(cmd):
    log.debug(cmd)
    if log.root.level <= 10:  # only show stdout/stderr if level is log.DEBUG
        redirect = None
    else:
        redirect = subprocess.DEVNULL
    subprocess.run(cmd, shell=True, stdout=redirect, stderr=redirect)


def rename_duplicates_in_list(input_list):
    """
    used to rename a dataframe's columns when multiple exist
    case insensitive
    """
    output_list = []
    elt_counts = {}
    for e in input_list:
        e_lower = e.lower()
        if e_lower not in elt_counts.keys():
            output_list.append(e)
            elt_counts[e_lower] = 1
        else:
            elt_counts[e_lower] += 1
            output_list.append(e + "_" + str(elt_counts[e_lower]))
    return output_list


def is_helper_func(arg):
    if isinstance(arg, list):
        if arg[0] == "HELPER_FUNCTION":
            return True
        else:
            raise ValueError(
                "This config file value should be a String or a HELPER_FUNCTION pattern:" + arg
            )
    return False


@lru_cache
def get_genome(fasta_path):
    return Fasta(fasta_path)


@lru_cache
def varank_to_vcf_coords(coord_conversion_file):
    """
    outside of helper class to avoid caching issues
    """
    id_to_coords = {}
    with open(coord_conversion_file, "r") as f:
        next(f)
        for l in f:
            l = l.strip().split("\t")
            id_to_coords[l[0]] = {
                "#CHROM": "chr" + l[1],
                "POS": l[2],
                "REF": l[3],
                "ALT": l[4],
            }
    return id_to_coords


def clean_string(s):
    """
    replace characters that will crash bcftools and/or cutevariant
    those in particular come from Varank files

    NB: the "fmt: off/on" comments are used to prevent black
    from making the replace dict into a one line mess
    """
    # fmt: off
    replace = {
        ";": ",",
        "“": '"',
        "”": '"',
        "‘": "'",
        "’": "'"
    }
    # fmt: on
    for k, v in replace.items():
        s = s.replace(k, v)
    return s


def create_vcf_header(
    input_path, config, sample_list, breakpoints=False, supplemental_header=[], fileformat="VCFv4.3"
):
    header = []
    header.append(f"##fileformat={fileformat}")
    header.append(f"##fileDate={time.strftime('%d/%m/%Y')}")
    header.append(f"##inputFile={os.path.abspath(input_path)}")
    header.append(f"##source=variantconvert from {config['GENERAL']['origin']}")
    header.append(f"##variantconvertVersion={variantconvert.__version__}")
    if "mode" in config["GENERAL"]:
        header.append(f"##variantconvertMode={config['GENERAL']['mode']}")

    # TODO: FILTER is not present in any tool implemented yet
    # so all variants are set to PASS
    if config["VCF_COLUMNS"]["FILTER"] != "" and config["GENERAL"]["origin"] != "AnnotSV":
        raise ValueError(
            "Filters are not implemented yet. "
            'Leave config["COLUMNS_DESCRIPTION"]["FILTER"] empty '
            "or whip the developer until he does it."
            "If you are trying to convert an annotSV file, "
            'use "annotsv" in the input file format argument'
        )
    header.append('##FILTER=<ID=PASS,Description="Passed filter">')

    if "ALT" in config["COLUMNS_DESCRIPTION"]:
        for key, desc in config["COLUMNS_DESCRIPTION"]["ALT"].items():
            header.append("##ALT=<ID=" + key + ',Description="' + desc + '">')

    if "INFO" in config["COLUMNS_DESCRIPTION"]:
        info_dic = config["COLUMNS_DESCRIPTION"]["INFO"]
        if breakpoints:
            if "SVTYPE" not in info_dic.keys():
                info_dic["SVTYPE"] = {
                    "Type": "String",
                    "Number": "1",
                    "Description": "Type of structural variant",
                }
            if "MATEDID" not in info_dic.keys():
                info_dic["MATEID"] = {
                    "Type": "String",
                    "Number": "1",
                    "Description": "ID of mate breakends",
                }

        for key, dic in info_dic.items():
            header.append(
                "##INFO=<ID="
                + key
                + ",Number="
                + dic["Number"]
                + ",Type="
                + dic["Type"]
                + ',Description="'
                + dic["Description"]
                + '">'
            )

    if supplemental_header != []:
        header += supplemental_header

    if "FORMAT" in config["COLUMNS_DESCRIPTION"]:
        for key, dic in config["COLUMNS_DESCRIPTION"]["FORMAT"].items():
            header.append(
                "##FORMAT=<ID="
                + key
                + ",Number=1,Type="
                + dic["Type"]
                + ',Description="'
                + dic["Description"]
                + '">'
            )

    header += config["GENOME"]["vcf_header"]
    header.append(
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + sample_list
        )
    )
    return header


def remove_decimal_or_strip(value):
    # To prevent bad format conversions by pandas with INFO field
    if is_float(value):
        value = str(float(value))
    elif is_int(value):
        value = str(value)
    if value.endswith(".0"):
        value = value[:-2]
    # forbidden to have spaces in INFO fields in VCF v4.2 (for IGV compatibility)
    value = value.replace(" ", "_")
    return value


def info_string_to_dict(info):
    """
    >>> info_string_to_dict("SVTYPE=duplication;SVLEN=35918771;END=73775259;SOMATIC;")
    {"SVTYPE":"duplication, "SVLEN": "35918771", "END": "73775259", "SOMATIC": None}
    """
    if ";" in info:
        data = info.split(";")
    else:
        data = [info]
    # if the INFO field mistakenly finished by a ';' then remove the final empty value from data
    if data[-1] == "":
        data = data[:-1]
    res = {}
    for pair in data:
        if pair.count("=") == 1:
            pair = pair.split("=")
            res[pair[0]] = pair[1]
        elif pair.count("=") == 0:
            res[pair] = None
        else:
            raise ValueError(
                f"info_string_to_dict(info) expects info to have the format 'key1=value1;key2=value2'. Got instead: {info}"
            )

    return res


def is_int(element: any) -> bool:
    if element == None:
        return False
    try:
        int(element)
        return True
    except ValueError:
        return False


def is_float(element: any) -> bool:
    if element == None:
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False


def all_equal(l: list) -> bool:
    """
    >>> all_equal(["a", "a", "a"])
    True
    >>> all_equal(["a", "b", "a"])
    False
    >>> all_equal([])
    True
    """
    try:
        if l.count(l[0]) == len(l):
            return True
    except IndexError:
        return True  # faster than checking every time if len(v) > 0
    return False
