# -*- coding: utf-8 -*-

import logging as log
import os
import pandas as pd
import re
import sys
import time
from natsort import index_natsorted

from variantconvert.converters.abstract_converter import AbstractConverter
from variantconvert.converters.vcf_from_annotsv_tools import check_config, merge_full_and_split

from variantconvert.commons import (
    create_vcf_header,
    info_string_to_dict,
    is_helper_func,
    is_int,
    is_float,
    remove_decimal_or_strip,
)

from variantconvert.helpers.annotsv_helper import AnnotSvHelper


class VcfFromAnnotsv(AbstractConverter):
    """
    Specificities compared to TSV:
    - vcf-like INFO field  in addition to other annotation columns
    - vcf-like FORMAT and <sample> fields
    - full/split annotations. Each variant can have one "full" and
    zero to many "split" annotations which result in additional lines in the file
    Very hard to deal with this with generic code --> it gets its own converter

    Sept 2022 update: add support for AnnotSV files obtained from bed files
    Those do not have REF, ALT, FORMAT and <sample_name> columns
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.config = check_config(self.config)

    def _build_input_dataframe(self):
        df = pd.read_csv(
            self.filepath,
            skiprows=self.config["GENERAL"]["skip_rows"],
            sep="\t",
            low_memory=False,
        )
        df.sort_values(
            [self.config["VCF_COLUMNS"]["#CHROM"], self.config["VCF_COLUMNS"]["INFO"]["SV_start"]],
            inplace=True,
        )
        df.reset_index(drop=True, inplace=True)

        sample_col = self.config["VCF_COLUMNS"]["SAMPLE"]
        if isinstance(sample_col, str) and sample_col != "":
            # avoid replacing "NA" sample by a dot
            df.loc[:, df.columns != sample_col] = df.loc[:, df.columns != sample_col].fillna(".")
            df.loc[:, sample_col] = df.loc[:, sample_col].fillna("NA")
        else:
            df.fillna(".", inplace=True)  # default empty value in VCF

        df = df.astype(str)
        log.debug(df)
        return df

    def _get_sample_list(self):
        if self.config["VCF_COLUMNS"]["FORMAT"] == "FORMAT":
            return self._get_sample_list_with_vcf_input()
        else:
            return self._get_sample_list_with_bed_input()

    def _get_sample_list_with_bed_input(self):
        samples_col = self.input_df[self.config["VCF_COLUMNS"]["SAMPLE"]]
        sample_list = []
        for cell in samples_col:
            if "," in cell:
                for sample in cell.split(","):
                    sample_list.append(sample)
            else:
                sample_list.append(cell)
        sample_list = list(set(sample_list))
        sample_list.sort()  # ensures output is always the same, despite using a set() above
        return sample_list

    def _get_sample_list_with_vcf_input(self):
        sample_list = []
        try:
            format_index = self.input_df.columns.to_list().index("FORMAT")
        except ValueError:
            raise ValueError(
                "ERROR: No FORMAT column detected. Variantconvert assumes you used a bed file as input for AnnotSV in that case. \nYou should use the 'annotsv3_from_bed.json' config, not 'annotsv3_from_vcf.json'"
            )

        for i in range(format_index + 1, self.input_df.columns.size):
            potential_sample_col = self.input_df[self.input_df.columns[i]]
            first_valid_index = potential_sample_col.first_valid_index()
            if self.is_sample_column(
                self.input_df.iloc[first_valid_index, i],
                self.input_df.iloc[first_valid_index, format_index],
            ):
                sample_list.append(self.input_df.columns[i])
            else:
                break

        log.debug(f"sample_list:{sample_list}")
        if not set(sample_list).issubset(self.input_df.columns):
            raise ValueError(
                f"When using an AnnotSV file generated from a VCF, all samples in {self.config['VCF_COLUMNS']['SAMPLE']} column are expected to have their own column in the input AnnotSV file"
            )
        return sample_list

    @staticmethod
    def is_sample_column(value, format):
        """
        Unperfect identification of sample columns. Assumed to be good enough within the context of AnnotSV files. Made to be used on columns directly after the FORMAT columns, knowing after sample columns should come the "Annotation" column (=full/split designation).
        Sample columns are arbitrarily identified because the "Samples" column
        """
        count = format.count(":")
        if count == 0:
            if format == "GT":
                # (.([|/].)+)  --> match if value is for example 0/1 or 1/0 or 1/2 or 0|1 or 0/1/2 etc...
                res = re.match(r"(.([|/].)+)", value)
                if res:
                    return True
                return False
            else:
                raise RuntimeError(
                    "Unable to determine which columns are sample columns.\nPlease verify that your input files are correct, then contact the developers at https://github.com/SamuelNicaise/variantconvert/issues"
                )
        else:
            if value.count(":") == count:
                return True
            return False

    def _build_input_annot_df(self) -> pd.DataFrame:
        """
        remove vcf base cols, FORMAT, Samples_ID, and each <sample> column
        expand 'INFO' field from origin vcf into new annotation columns
        perform some checks
        """
        columns_to_drop = [v for v in self.sample_list]
        columns_to_drop += [v for v in self.main_vcf_cols]
        columns_to_drop.append(self.config["VCF_COLUMNS"]["SAMPLE"])
        columns_to_drop.append(self.config["VCF_COLUMNS"]["FORMAT"])

        # TODO: verify all booleans on config load
        if self.config["GENERAL"].get("keep_info", False) in (True, "true", "True"):
            self.keep_info = True
            if self.config["GENERAL"]["mode"] == "combined":
                self.original_info_col = self.input_df[
                    [
                        self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"],
                        self.config["VCF_COLUMNS"]["INFO"]["INFO"],
                    ]
                ]
                self.original_info_col.set_index(
                    self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"], inplace=True
                )
                self.original_info_col = self.original_info_col["INFO"].to_dict()
            else:
                self.original_info_col = self.input_df[self.config["VCF_COLUMNS"]["INFO"]["INFO"]]
                self.original_info_col = self.original_info_col.to_dict()
        else:
            self.keep_info = False
        # columns_to_drop.append(self.config["VCF_COLUMNS"]["INFO"]["INFO"])

        df = self.input_df.copy()
        for col in columns_to_drop:
            try:
                df = df.drop([col], axis=1)
            except KeyError:
                log.debug(f"Failed to drop column: {col}")

        df = df.replace(";", "|", regex=True)  # any ';' in annots will ruin the vcf INFO field

        # TODO: check if CHROM col is in compliance with config ref genome (chrX or X)
        # if self.config["GENOME"]["vcf_header"][0].startswith("##contig=<ID=chr"):
        #     if not chrom.startswith

        if (
            self.config["VCF_COLUMNS"]["INFO"]["SVTYPE"] == ""
            or self.config["VCF_COLUMNS"]["INFO"]["SVTYPE"] not in df.columns
        ):
            raise ValueError(
                "SV_type column is required to turn an AnnotSV file into a VCF. Check if SV_type col is set in config or missing in your file.\n"
                + "If you generated your AnnotSV file from a bed, AnnotSV option -svtBEDcol is required."
            )
        return df

    def _build_info_dic(self):
        """
        Output: dictionary with key: annotsv_ID ; value: a key-value dictionary of all annotations
        This will be used to write the INFO field
        """
        self.input_annot_df = self._build_input_annot_df()
        annots_dic = {}
        supplemental_info_fields = []

        # method of versions < 2.0.0: merge full and split annotations into one variant
        if self.config["GENERAL"]["mode"] == "combined":
            id_col = self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"]
            for variant_id, df_variant in self.input_annot_df.groupby(id_col):
                merged_annots = merge_full_and_split(self.config, self.sample_list, df_variant)
                annots_dic[variant_id] = merged_annots

                if self.keep_info:
                    info_dict = info_string_to_dict(self.original_info_col[variant_id])
                    for annot in info_dict:
                        if annot not in annots_dic[variant_id]:
                            annots_dic[variant_id][annot] = info_dict[annot]
                            if annot not in supplemental_info_fields:
                                supplemental_info_fields.append(annot)

        else:
            for variant_idx, df_variant in self.input_annot_df.iterrows():
                annots_dic[variant_idx] = df_variant.to_dict()

                if variant_idx == 0:
                    log.debug(f"annots dic first variant before adding INFO: {annots_dic[0]}")

                if self.keep_info:
                    info_dict = info_string_to_dict(self.original_info_col[variant_idx])
                    for annot in info_dict:
                        if annot not in annots_dic[variant_idx]:
                            annots_dic[variant_idx][annot] = info_dict[annot]
                            if annot not in supplemental_info_fields:
                                supplemental_info_fields.append(annot)
                else:
                    annots_dic[variant_idx].pop("INFO", None)

        if self.keep_info:
            self.supplemental_header = self._build_supplemental_header(
                annots_dic, supplemental_info_fields
            )
        else:
            self.supplemental_header = self._build_supplemental_header(annots_dic, [])

        # log.debug(f"annots dic first variant after adding info: {annots_dic[0]}")
        return annots_dic

    def _build_supplemental_header(self, annots_dic, additional_info_fields):
        """
        Create header strings for each value in additional_info_fields
        Types are inferred from the values in annots_dic

        Args:
            annots_dic (dict): has the following structure:
                                {
                                    variant1: {key1: val1, key2: val2},
                                    variant2: {key1: val1, key2: val2}
                                }
            additional_info_fields (set[str]): new info fields who need to be added to header.
                                            This list will be automatically completed with columns that
                                            are in annots_dic but not defined in config["COLUMNS_DESCRIPTION]

        Returns:
            list[str]: header strings, one for each additional_info_fields.
            Doesn't include newlines, as expected by commons.create_vcf_header(args)
        """
        supplemental_header = []
        known_descriptions = set(self.config["COLUMNS_DESCRIPTION"]["INFO"].keys())
        missing_annots = []

        # this function was originally made to rebuild the lost header of VCF annotations in "INFO" column in VCF>AnnotSV>VCF conversions.
        # at this point, we check if all INFO fields were defined in config.
        # if not, add a default header for them too.
        for variant, dic in annots_dic.items():
            for annot in dic:
                if (
                    annot not in known_descriptions
                    and annot not in missing_annots
                    and annot not in ["SV_end", "SV_length", "SV_type"]
                ):
                    missing_annots.append(annot)

        missing_annots = additional_info_fields + [
            v for v in missing_annots if v not in additional_info_fields
        ]
        # filter a second time because additional info fields may be have duplicates with known descriptions
        missing_annots = [v for v in missing_annots if v not in known_descriptions]
        log.debug(f"known desc:{missing_annots}")

        for field in missing_annots:
            # infer type
            if all([is_int(v.get(field, None)) for v in annots_dic.values()]):
                field_type = "Integer"
            elif all([is_float(v.get(field, None)) for v in annots_dic.values()]):
                field_type = "Float"
            else:
                field_type = "String"

            # infer number
            for info in annots_dic.values():
                if field in info:
                    if info[field] == None:
                        number = 0
                        field_type = "Flag"
                    else:
                        # used to infer: number = info[field].count(",") + 1
                        # instead: always put "." as it can change for each variant depending on the number of split annots
                        number = "."

            if field in additional_info_fields:
                description = (
                    "Imported from the INFO field of the original VCF before AnnotSV annotation"
                )
            else:
                description = "Imported from AnnotSV"

            supplemental_header.append(
                "##INFO=<ID="
                + field
                + ",Number="
                + str(number)
                + ",Type="
                + field_type
                + ',Description="'
                + description
                + '">'
            )
        return supplemental_header

    def _get_main_vcf_cols(self):
        """
        Some columns (cols_to_init_now) are not directly linked to a Varank TSV column.
        They can be filled later, usually with a HELPER FUNC, but still have to be initiated in the data frame.
        """
        main_cols = []
        cols_to_init_now = []

        for vcf_col in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]:
            tsv_col = self.config["VCF_COLUMNS"][vcf_col]
            if not isinstance(tsv_col, str):
                cols_to_init_now.append(vcf_col)
            elif tsv_col == "":
                cols_to_init_now.append(vcf_col)
            else:
                vcf_col = tsv_col
            main_cols.append(vcf_col)

        # adding all missing columns at once to avoid a PerformanceWarning
        self.input_df = pd.concat([self.input_df, pd.DataFrame(columns=cols_to_init_now)])
        self.input_df[cols_to_init_now] = ["." for i in range(len(cols_to_init_now))]

        return main_cols

    def _get_variant_info_dict(self, input_dic, helper):
        """
        input_dic is the info_dict entry corresponding to the current line being written
        It has already some modifications done such as transforming the "INFO" column into separate annotations

        Returns the final annotations dict that will be written into the output vcf for the current line
        """
        # not info fields
        ignored_cols = self.main_vcf_cols + self.sample_list
        # END, SVLEN and SVTYPE are reserved INFO keywords for SVs per VCF 4.3 specification
        # They are renamed by creating new columns through config, and deleting the old ones
        ignored_cols += ["SV_end", "SV_length", "SV_type"]
        for config_key, config_val in self.config["VCF_COLUMNS"]["INFO"].items():
            if is_helper_func(config_val):
                func = helper.get(config_val[1])
                args = [input_dic.get(c, None) for c in config_val[2:]]
                result = func(*args)
                input_dic[config_key] = result

        cleaned_info_dic = {}
        for k, v in input_dic.items():
            if k not in ignored_cols:
                if (
                    self.config["GENERAL"]["keep_empty_info"] in (False, "false", "False")
                    and v == "."
                ):
                    continue
                if v != None:
                    if isinstance(v, str):
                        v = v.replace(";", "|")
                    cleaned_info_dic[k] = remove_decimal_or_strip(v)
                else:
                    cleaned_info_dic[k] = None
        return cleaned_info_dic

    def convert(self, tsv, output_path):
        """
        Creates and fill the output file.

        For each annotSV_ID ; fetch all related lines of annotations in key value dics.
        A function identifies and merges the annotations.
        then we build a dictionary with 1 dictionary per annotSV_ID containing all the annotations
        This is then used to make the header and fill the INFO field.
        """
        log.debug("Converting to vcf from tsv using config: " + self.config_filepath)

        self.filepath = tsv
        helper = AnnotSvHelper(self.config)

        self.input_df = self._build_input_dataframe()
        self.sample_list = self._get_sample_list()
        self.main_vcf_cols = self._get_main_vcf_cols()

        info_input_dic = self._build_info_dic()

        # create the vcf
        with open(output_path, "w") as vcf:
            vcf_header = create_vcf_header(
                tsv,
                self.config,
                self.sample_list,
                supplemental_header=self.supplemental_header,
                fileformat="VCFv4.2",  # because no spaces in INFO fields for IGV support
            )
            for l in vcf_header:
                vcf.write(l + "\n")

            # id_col = self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"]
            self.input_df = self.input_df.iloc[
                index_natsorted(self.input_df[self.config["VCF_COLUMNS"]["#CHROM"]])
            ]

            already_seen_variants = set()  # for combined mode

            for row_idx, row_data in self.input_df.iterrows():
                df_variant = row_data.to_dict()  # for ease of dev
                if self.config["GENERAL"]["mode"] == "combined":
                    variant_idx = df_variant[self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"]]
                    if variant_idx in already_seen_variants:
                        continue
                    else:
                        already_seen_variants.add(variant_idx)
                else:
                    variant_idx = row_idx

                if (
                    self.config["GENERAL"]["mode"] == "full"
                    and df_variant[self.config["VCF_COLUMNS"]["INFO"]["Annotation_mode"]] != "full"
                ):
                    continue

                # fill columns that need a helper func
                for config_key, config_val in self.config["VCF_COLUMNS"].items():
                    if config_key == "INFO":
                        continue  # dealt with by _get_variant_info_dict()
                    elif config_key == "FILTER" and config_val == "":
                        df_variant[config_key] = "PASS"
                    elif is_helper_func(config_val):
                        func = helper.get(config_val[1])
                        args = [df_variant.get(c, None) for c in config_val[2:]]
                        result = func(*args)
                        df_variant[config_key] = result

                main_cols = "\t".join([df_variant[c] for c in self.main_vcf_cols])
                vcf.write(main_cols + "\t")

                # INFO field
                info_list = []
                for k, v in self._get_variant_info_dict(
                    info_input_dic[variant_idx], helper
                ).items():
                    if v != None:
                        info_list.append(k + "=" + v)
                    else:
                        info_list.append(k)  # deal with INFO flags
                vcf.write(";".join(info_list) + "\t")

                # FORMAT and samples
                if self.config["VCF_COLUMNS"]["FORMAT"] != "":
                    sample_cols = "\t".join(
                        [df_variant[self.config["VCF_COLUMNS"]["FORMAT"]]]
                        + [df_variant[c] for c in self.sample_list]
                    )
                else:
                    sample_cols = "GT"
                    samples_with_variant = df_variant[self.config["VCF_COLUMNS"]["SAMPLE"]].split(
                        ","
                    )
                    for sample in self.sample_list:
                        if sample in samples_with_variant:
                            sample_cols += f"\t{self.config['GENERAL']['default_present_genotype']}"
                        else:
                            sample_cols += f"\t{self.config['GENERAL']['default_absent_genotype']}"
                vcf.write(sample_cols)
                vcf.write("\n")
