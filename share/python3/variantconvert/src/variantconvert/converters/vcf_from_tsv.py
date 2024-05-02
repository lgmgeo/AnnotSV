# -*- coding: utf-8 -*-

import logging as log
import os
import pandas as pd
import sys

from converters.abstract_converter import AbstractConverter

sys.path.append("..")
from commons import create_vcf_header, is_helper_func, clean_string
from variantconvert.helpers.helper_functions import HelperFunctions


class VcfFromTsv(AbstractConverter):
    def _init_dataframe(self):
        self.df = pd.read_csv(
            self.filepath,
            skiprows=self.config["GENERAL"]["skip_rows"],
            sep="\t",
            low_memory=False,
        )
        self.df.sort_values(
            [self.config["VCF_COLUMNS"]["#CHROM"], self.config["VCF_COLUMNS"]["POS"]],
            inplace=True,
        )
        self.df.reset_index(drop=True, inplace=True)
        self.df.fillna(".", inplace=True)
        log.debug(self.df)
        self.df[self.UNIQUE_ID] = self.df.apply(
            lambda row: self._get_unique_variant_id(row), axis=1
        )
        if self.config["VCF_COLUMNS"]["SAMPLE"] != "":
            self.df[self.config["VCF_COLUMNS"]["SAMPLE"]] = self.df.apply(
                lambda row: self._bwamem_name_bugfix(row), axis=1
            )
        log.debug(self.df)

    def _get_sample_list(self):
        # is the file multisample?
        if self.config["VCF_COLUMNS"]["SAMPLE"] != "":
            sample_list = []
            for sample in self.df[self.config["VCF_COLUMNS"]["SAMPLE"]].unique():
                sample_list.append(sample)
            return sample_list
        else:
            return [os.path.basename(self.output_path)]

    def _bwamem_name_bugfix(self, row):
        """remove .bwamem from the end of sample names if needed"""
        name = row[self.config["VCF_COLUMNS"]["SAMPLE"]]
        if name.endswith(".bwamem") and name != ".bwamem":
            return name[0:-7]
        else:
            return name

    def _get_unique_variant_id(self, row):
        var_id = []
        for col in self.config["GENERAL"]["unique_variant_id"]:
            var_id.append(str(row[col]))
        return "_".join(var_id)

    def _get_unique_id_to_index_list(self, data):
        id_dic = {}
        for k, v in data[self.UNIQUE_ID].items():
            if v not in id_dic:
                id_dic[v] = [k]
            else:
                id_dic[v].append(k)
        return id_dic

    def convert(self, tsv, output_path):
        log.debug("Converting to vcf from annotSV using config: " + self.config_filepath)

        self.UNIQUE_ID = "__!UNIQUE_VARIANT_ID!__"
        self.filepath = tsv
        self.output_path = output_path
        self._init_dataframe()
        sample_list = self._get_sample_list()
        log.debug(f"sample_list:{sample_list}")
        helper = HelperFunctions(self.config)

        with open(output_path, "w") as vcf:
            vcf_header = create_vcf_header(tsv, self.config, sample_list)
            for l in vcf_header:
                vcf.write(l + "\n")

            data = self.df.astype(str).to_dict()
            # In Decon (and maybe others), TSV are given as a list of variant-sample associations
            # so the same variant can be on multiple TSV lines
            # __!UNIQUE_VARIANT_ID!__ allows to identify variants and only add them to the VCF once
            already_seen_variants = set()
            unique_id_to_index_list = self._get_unique_id_to_index_list(data)
            for i in range(len(data[self.config["VCF_COLUMNS"]["#CHROM"]])):
                if data[self.UNIQUE_ID][i] in already_seen_variants:
                    continue

                line = ""

                for vcf_col in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL"]:
                    col = self.config["VCF_COLUMNS"][vcf_col]
                    if is_helper_func(col):
                        # col[1] is a function name, col[2] its list of args
                        # the function named in col[1] has to be callable from this module
                        func = helper.get(col[1])
                        args = [data[c][i] for c in col[2:]]
                        line += func(*args) + "\t"
                    elif col == "":
                        line += ".\t"
                    else:
                        line += data[col][i] + "\t"

                # Cutting-edge FILTER implementation
                line += "PASS\t"

                info_field = []
                for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["INFO"].items():
                    if is_helper_func(tsv_col):
                        func = helper.get(tsv_col[1])
                        args = [data[c][i] for c in tsv_col[2:]]
                        s = vcf_col + "=" + func(*args)
                    else:
                        s = vcf_col + "=" + data[tsv_col][i]
                    info_field.append(clean_string(s))
                line += ";".join(info_field) + "\t"

                vcf_format_fields = []
                tsv_format_fields = []
                for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["FORMAT"].items():
                    vcf_format_fields.append(vcf_col)
                    tsv_format_fields.append(tsv_col)
                line += ":".join(vcf_format_fields) + "\t"

                # monosample input
                if len(sample_list) == 1:
                    sample_field = []
                    for index in unique_id_to_index_list[data[self.UNIQUE_ID][i]]:
                        for key, val in self.config["VCF_COLUMNS"]["FORMAT"].items():
                            if key == "GT" and val == "":
                                sample_field.append("0/1")
                                continue
                            sample_field.append(data[val][index])
                        line += ":".join(sample_field)
                    # TODO: deal with monosample files with variants that are not associated to any sample

                # multisample input
                else:
                    sample_field_dic = {}
                    # If the variant exists in other lines in the source file, fetch their sample data now
                    for index in unique_id_to_index_list[data[self.UNIQUE_ID][i]]:
                        sample_field = []
                        for key, val in self.config["VCF_COLUMNS"]["FORMAT"].items():
                            if key == "GT" and val == "":
                                sample_field.append("0/1")
                                continue
                            sample_field.append(data[val][index])
                        sample_field_dic[data[self.config["VCF_COLUMNS"]["SAMPLE"]][index]] = (
                            ":".join(sample_field)
                        )

                    for sample in sample_list:
                        if sample in sample_field_dic:
                            line += sample_field_dic[sample] + "\t"
                        else:
                            if "GT" in self.config["VCF_COLUMNS"]["FORMAT"]:
                                if len(self.config["VCF_COLUMNS"]["FORMAT"]) == 1:
                                    # there's only GT. Avoid adding a trailing ":"
                                    empty = "./."
                                else:
                                    empty = "./.:" + ":".join(
                                        [
                                            "."
                                            for i in range(
                                                len(self.config["VCF_COLUMNS"]["FORMAT"]) - 1
                                            )
                                        ]
                                    )
                            else:
                                empty = ":".join(
                                    [
                                        "."
                                        for i in range(
                                            len(self.config["VCF_COLUMNS"]["FORMAT"]) - 1
                                        )
                                    ]
                                )
                            line += empty + "\t"
                    already_seen_variants.add(data[self.UNIQUE_ID][i])
                    line = line.rstrip("\t")

                vcf.write(line + "\n")
