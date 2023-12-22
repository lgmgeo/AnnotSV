# -*- coding: utf-8 -*-
"""
"""

import filecmp
import json
import logging as log
import os
import pytest
import sys

from os.path import join as osj

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from variantconvert.config import change_config, get_nested_dic, update_nested_dic


FAKE_CONFIG = {
    "GENERAL": {
        "origin": "AnnotSV",
        "skip_rows": 0,
    },
    "GENOME": {
        "assembly": "GRCh37",
        "path": "human_g1k_v37.fasta",
        "vcf_header": [
            "##contig=<ID=1,length=249250621,assembly=GRCh37>",
            "##contig=<ID=2,length=243199373,assembly=GRCh37>",
        ],
    },
}

NEW_VARS_1 = {"GENOME": {"assembly": "hg19"}}
CONTROL_1 = {
    "GENERAL": {
        "origin": "AnnotSV",
        "skip_rows": 0,
    },
    "GENOME": {
        "assembly": "hg19",
        "path": "human_g1k_v37.fasta",
        "vcf_header": [
            "##contig=<ID=1,length=249250621,assembly=GRCh37>",
            "##contig=<ID=2,length=243199373,assembly=GRCh37>",
        ],
    },
}

NEW_VARS_2 = {
    "GENERAL": {"origin": "tests", "skip_rows": "no longer an int", "new_key": "new_val"},
    "GENOME": {"assembly": "hg19"},
    "NEW_SECTION": "other_val",
}
CONTROL_2 = {
    "GENERAL": {
        "origin": "tests",
        "skip_rows": "no longer an int",
        "new_key": "new_val",
    },
    "GENOME": {
        "assembly": "hg19",
        "path": "human_g1k_v37.fasta",
        "vcf_header": [
            "##contig=<ID=1,length=249250621,assembly=GRCh37>",
            "##contig=<ID=2,length=243199373,assembly=GRCh37>",
        ],
    },
    "NEW_SECTION": "other_val",
}


def test_get_nested_dic():
    dic = {}
    target_key = "COLUMNS.INFO.ID.Test"
    expected = {"COLUMNS": {"INFO": {"ID": {"Test": ""}}}}
    res1 = get_nested_dic(dic, target_key)
    assert res1 == expected

    target_key = "COLUMNS.FILTER"
    default_val = "Pass"
    # fmt: off
    expected = {'COLUMNS': {
                        'INFO': {
                            'ID': {
                                'Test': ''
                            }
                        },
                        'FILTER': 'Pass'
                    } 
                }
    # fmt: on
    res2 = get_nested_dic(res1, target_key, default_value=default_val)
    assert res2 == expected


def test_update_nested_dic():
    assert CONTROL_1 == update_nested_dic(FAKE_CONFIG, NEW_VARS_1)
    assert CONTROL_2 == update_nested_dic(FAKE_CONFIG, NEW_VARS_2)


def test_change_config(tmp_path):

    new_vars = [NEW_VARS_1, NEW_VARS_2]
    controls = [CONTROL_1, CONTROL_2]

    initial_json = osj(tmp_path, "initial.json")
    test_json = osj(tmp_path, "test.json")

    for new, control in zip(new_vars, controls):
        # write two identical initial files
        with open(initial_json, "w") as f:
            f.write(json.dumps(FAKE_CONFIG, indent="\t"))
        with open(test_json, "w") as f:
            f.write(json.dumps(FAKE_CONFIG, indent="\t"))

        # update test.json, then compare with initial.json
        change_config(test_json, new)
        assert filecmp.cmp(initial_json, test_json)
