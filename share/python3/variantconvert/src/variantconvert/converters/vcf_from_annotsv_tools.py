# -*- coding: utf-8 -*-

"""
Splitting annotsv converter as it keeps getting bigger.
This file is dedicated to functions for type 3) AnnotSV > VCF conversion. See below for explanation about types 1,2,3.


With this given AnnotSV input
AnnotSV_ID	SV_chrom	SV_start	SV_end	SV_length	Variant_type	Annotation_mode	Gene_name	DDD_HI_percent
10_46976157_47590995__1	10	46976157	47590995		DUP	full	AGAP9;ANTXRLP1;ANXA8	91.07
10_46976157_47590995__1	10	46976157	47590995		DUP	split	AGAP9	88.1
10_46976157_47590995__1	10	46976157	47590995		DUP	split	ANTXRLP1	
10_46976157_47590995__1	10	46976157	47590995		DUP	split	ANXA8	76.24


AnnotSV converter has 3 possible outputs
1) keep only full

##DDD_HI_percent type=Float
#CHROM	POS	REF	ALT	INFO
chr10	46976157	N	<DUP>	Annotation=full ; DDD_HI_percent=91.07

2) keep full + one line per split annotation

##DDD_HI_percent type=Float
#CHROM	POS	REF	ALT	INFO
chr10	46976157	N	<DUP>	Annotation=full ; DDD_HI_percent=91.07
chr10	46976157	N	<DUP>	Annotation=split ; DDD_HI_percent=88.1
chr10	46976157	N	<DUP>	Annotation=split
chr10	46976157	N	<DUP>	Annotation=split ; DDD_HI_percent=76.24

3) variantconvert version < 2.0 : 1 line per variant with all annotations packed together
##DDD_HI_percent type=String
#CHROM	POS	REF	ALT	INFO
chr10	46976157	N	<DUP>	Annotation=full|split|split|split; DDD_HI_percent=91.07|88.1|.|76.24

Two changes for type 3) in variantconvert >=2.0 : 
- replace the variant convert "|" with ","
- replace the "," in annotations with "|"
"""

import logging as log
import pandas as pd

from variantconvert.commons import all_equal, remove_decimal_or_strip


def check_config(config: dict) -> dict:
    ACCEPTED_MODES = ("full", "full&split", "combined")
    if config["GENERAL"]["mode"] not in ACCEPTED_MODES:
        raise ValueError(
            f"Unexpected value in json config['GENERAL']['mode']: {config['GENERAL']['mode']} -- accepted values: {ACCEPTED_MODES}"
        )
    # NB: returns config instead of void as some inferred changes may be automatically made in the future
    return config


def merge_full_and_split(config: dict, sample_list: list[str], df: pd.DataFrame, sep=",") -> dict:
    """
    input: df of a single annotSV_ID ; containing only annotations (no sample/FORMAT data)
    it can contain full and/or split annotations

    returns a single line dataframe with all annotations merged properly
    """
    # Do not keep 'base vcf col' in info field
    df = df.loc[
        :,
        [
            cols
            for cols in df.columns
            if cols not in ["ID", "REF", "ALT", "QUAL", "FILTER"] + sample_list
        ],
    ]
    annots = {}
    dfs = {}
    for typemode, df_type in df.groupby(config["VCF_COLUMNS"]["INFO"]["Annotation_mode"]):
        if typemode not in ("full", "split"):
            raise ValueError("Annotation type is assumed to be only 'full' or 'split'")
        dfs[typemode] = df_type

    # deal with full
    if "full" not in dfs.keys():
        # still need to init columns
        if "split" not in dfs.keys():
            log.warning(
                "Input does not include AnnotSV's 'Annotation_mode' column. This is necessary to know how to deal with annotations. The INFO field will be empty."
            )
            return {}
        # do not initiate an
        # for ann in dfs["split"].columns:
        #     annots[ann] = "."
    else:
        if len(dfs["full"].index) > 1:
            raise ValueError(
                "Each variant is assumed to only have one single line of 'full' annotation"
            )
        # remove float decimal full rows
        for ann in dfs["full"].columns:
            annots[ann] = remove_decimal_or_strip(dfs["full"].loc[df.index[0], ann])

    # deal with split
    if "split" not in dfs.keys():
        return annots
    # each info field split
    for ann in dfs["split"].columns:
        transform = []
        # list of all values in each columns
        for splitval in dfs["split"][ann].tolist():
            transform.append(remove_decimal_or_strip(splitval))

        # don't report split infos if there are ALL undefined (=".")
        if all(eq == "." for eq in transform):
            continue
        else:
            try:
                values = [annots[ann]]
                values.extend(transform)
            except KeyError as err:
                # err.args[0] returns the missing key from KeyError exception
                log.debug(
                    f"Missing key {err.args[0]} - assuming variant has no 'full' annotation, only 'split'"
                )
                values = transform

        # merge annots who are all identical (SVLEN...)
        if all_equal(values):
            annots[ann] = values[0]
        else:
            annots[ann] = sep.join(values)

    return annots
