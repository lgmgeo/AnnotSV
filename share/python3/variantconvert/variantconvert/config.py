# -*- coding: utf-8 -*-
"""
Use example for testing
variantconvert config --set GENOME.assembly=swag GENOME.path=/home1/DB/STARK/genomes/current/hg19.fa GENERAL.origin=myapp
"""

import glob
import json
import logging as log
import os
import sys

from os.path import join as osj
from pyfaidx import Fasta

sys.path.append(os.path.join(os.path.dirname(__file__), "."))
from commons import set_log_level


def change_config(config_path, new_vars):
    with open(config_path, "r") as config:
        data = json.load(config)

    updated_data = update_nested_dic(data, new_vars)  # requires Python >= 3.5

    with open(config_path, "w") as config:
        config.write(json.dumps(updated_data, indent="\t"))


def update_nested_dic(data, new_vars):
    for k, v in new_vars.items():
        if isinstance(v, dict):
            data[k] = update_nested_dic(data.get(k, {}), v)
        else:
            data[k] = v
    return data


def get_nested_dic(dictionary, key_tree, sep=".", default_value=""):
    """
    >>> get_nested_dic({}, VCF_COLUMNS.INFO.ID)
    {"VCF_COLUMNS": {"INFO": {"ID": ""}}}

    Returns a nested dictionary containing the desired key_tree
    Value of the deepest key is initialized at default_value
    """
    keys = [key for key in key_tree.split(sep)]
    tmp_dic = dictionary
    for i in range(len(keys)):
        if keys[i] not in dictionary:
            if i == len(keys) - 1:
                tmp_dic[keys[i]] = default_value
            else:
                tmp_dic[keys[i]] = {}
                tmp_dic = tmp_dic[keys[i]]
        else:
            tmp_dic = tmp_dic[keys[i]]
    return dictionary


def fill_genome_header(config_path):
    """
    Using GENOME["assembly"] and GENOME["path"] from config_path, update GENOME["vcf_header"] in config_path.
    """
    with open(config_path, "r") as f:
        config = json.load(f)
    assembly = config["GENOME"]["assembly"]
    genome_path = config["GENOME"]["path"]
    fasta = Fasta(genome_path)

    vcf_header = []
    for contig in fasta:
        vcf_header.append(
            f"##contig=<ID={contig.long_name.split()[0]},length={len(contig)},assembly={assembly}>"
        )
    vcf_header.append(f"##reference=file://{genome_path}")
    log.debug(f"Built vcf_header:{vcf_header}")

    new_vars = {"GENOME": {"vcf_header": vcf_header}}
    change_config(config_path, new_vars)


def main_config(args):
    set_log_level(args.verbosity)

    if args.configFiles == "<script_dir>/configs/*":  # default value
        target_files = glob.glob(osj(os.path.dirname(__file__), "..", "configs", "*.json"))
    else:
        target_files = [v for v in args.configFiles]
    for target in target_files:
        if not os.path.exists(target):
            raise FileNotFoundError(target)
    log.debug(f"Config will be applied to files: {target_files}")

    new_vars = {}
    if args.set != None:
        for key_val in args.set:
            if key_val.count("=") != 1:
                raise RuntimeError(
                    "Use the following format: --set key1=value1 key2=value2 key3=value3"
                )
            items = key_val.split("=")
            key = items[0]
            val = items[1]
            new_vars = get_nested_dic(new_vars, key, default_value=val)

    log.debug(f"Config new_vars: {new_vars}")

    for config in target_files:
        if new_vars != {}:
            change_config(config, new_vars)

        if args.fill_genome_header == True:
            fill_genome_header(config)
