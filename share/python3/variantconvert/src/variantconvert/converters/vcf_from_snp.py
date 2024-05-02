"""
This class transforms raw data coming from SNP analysis into a vcf file:

exemple:

Raw_data input_file

Index   Name    Address Chr Position    ... Sample1.Top Alleles ... Sample2.Top Alleles ... n.Top Alleles
1       rs0001  3568    5   8759622     ... GG  ... GA  ... n.haplotypes
2       rs0002  4567    2   8972124     ... TA  ... CT  ... n.haplotypes
3       rs0045  7778    8   2655891     ... C-  ... CC  ... n.haplotypes

vcf output_file

##META_DATA
#CHROM  POS     ID      REF ALT QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 n.Samples
chr5    8759622 rs0001  G   A   .       PASS    .       GT      0/0     0/1     n.genotype
chr2    8972124 rs0002  T   A,C .       PASS    .       GT      0/1     1/0     n.genotype
chr8    2655891 rs0045  GC  G   .       PASS    .       GT      1/1     0/0     n.genotype

The vcf are unphased but you can use Beagle2.3 to phase your data

@Author: Elise Verin
"""

import logging as log
import pandas as pd
import os
import sys
import re

from converters.abstract_converter import AbstractConverter

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from variantconvert.helpers.helper_functions import HelperFunctions
from commons import create_vcf_header, is_helper_func
from variant import Variant


class VcfFromSnp(AbstractConverter):
    def _init_dataframe(self):
        self.snp_data = pd.read_csv(self.filepath, sep="\t", index_col=0)
        self.snp_data.reset_index(drop=True, inplace=True)

    def _get_sample_id(self):
        self.sample_list = []
        self.sample_alt = {}
        for i in range(self.snp_data.shape[1]):
            string = self.snp_data.columns[i]

            found_allele = re.search("Top Alleles", string)

            if found_allele != None:
                row_alt = {}
                # keep top allele from raw and split them
                split_allele = string.split(".")
                self.sample_list.append(split_allele[0])

                for a in range(self.snp_data.shape[0]):
                    # For each id find we keep all haplotype
                    row_alt[a] = self.snp_data.iloc[a, i]

                self.sample_alt[split_allele[0]] = row_alt

        return self.sample_list

    def manage_alt(self, row, ref):
        """
        The alt are generated based on ref and haplotype contained in row.
        When we detect a deletion we transform ref, all alleles and alt (if alt is not empty)

        Before deletion
        ref = A
        del_variant =T
        alt = T,C
        allele_1 = T
        allele_2 = C

        After we detect deletion
        ref= TA
        del_variant = T
        alt = TT,TC,T
        allele_1 = TT
        allele_2 = TC
        """

        alt = ""
        self.del_variant = False

        for i in self.sample_alt:
            values_sample = self.sample_alt.get(i)
            all_1 = values_sample.get(row)[0]
            all_2 = values_sample.get(row)[1]

            if all_1 == "-" or all_2 == "-":
                if self.del_variant == False:
                    self.del_variant = True
                    ref = self.nuc_del + ref

                    if alt != "":
                        alt_nw = ""
                        if alt == ".":
                            alt = ""
                        for str_alt in alt:
                            if str_alt != ",":
                                alt_nw = alt_nw + self.nuc_del + str_alt
                            else:
                                alt_nw = alt_nw + str_alt

                        alt = alt_nw

            if self.del_variant == True:
                if all_1 == "-":
                    all_1 = ""
                if all_2 == "-":
                    all_2 = ""
                all_1 = self.nuc_del + all_1
                all_2 = self.nuc_del + all_2

            alt = self.generate_alt(alt, all_1, all_2, ref)
        return alt

    def generate_alt(self, alt, all_1, all_2, ref):
        """For each allele we compare what ref and alt contain"""

        if alt == "":
            if all_1 == ref and all_2 == ref and alt.find(all_1) == -1 and alt.find(all_2) == -1:
                alt = "."

            if all_1 != ref and alt.find(all_1) == -1:
                alt = all_1

            if all_2 != ref and alt.find(all_2) == -1:
                alt = all_2

        elif alt != "":
            if all_1 != ref and alt.find(all_1) == -1:
                if alt.find(".") == 0:
                    alt = all_1
                else:
                    alt = alt + "," + all_1

            if all_2 != ref and alt.find(all_2) == -1:
                if alt.find(".") == 0:
                    alt = all_2
                else:
                    alt = alt + "," + all_2

        if self.del_variant == True and len(all_1) == 1 and alt != "":
            if not self.search_del(alt):
                alt = alt + "," + all_1

        elif self.del_variant == True and len(all_2) == 1 and alt != "":
            if not self.search_del(alt):
                alt = alt + "," + all_2

        return alt

    def search_del(self, alt):
        split_alt = alt.split(",")
        for i in split_alt:
            if len(i) == 1:
                return True
        return False

    def _define_gt(self, row, ref, sample):
        """If the nucleotide from the haplotype are found in REF we have a 0 and if it's found in alt it's 1. "/" correspond unphased genotype"""
        gt_samples = {}

        gt = ""
        value_sample = self.sample_alt[sample][row]

        for a in range(2):
            allele = str(value_sample[a])

            if allele == ref and gt == "":
                gt = "0"

            elif allele != ref and gt == "":
                gt = "1"

            elif allele == ref and gt != "":
                gt = gt + "/0"

            elif allele != ref and gt != "":
                gt = gt + "/1"

        gt_samples["GT"] = gt

        return gt_samples

    def convert(self, input_path, output_path):
        log.debug("Converting to vcf from SNP using config: " + self.config_filepath)

        self.filepath = input_path
        self.output_path = output_path
        self._init_dataframe()

        sample_list = self._get_sample_id()
        helper = HelperFunctions(self.config)

        with open(output_path, "w") as vcf:
            vcf_header = create_vcf_header(input_path, self.config, sample_list, True)

            for l in vcf_header:
                vcf.write(l + "\n")

            data = self.snp_data.astype(str).to_dict()

            error_reff = 0
            for i in range(self.snp_data.shape[0]):
                var = Variant()
                self.nuc_del = ""

                for vcf_col in ["#CHROM", "POS", "ID", "REF", "QUAL"]:
                    col = self.config["VCF_COLUMNS"][vcf_col]

                    if is_helper_func(col):
                        func = helper.get(col[1])
                        args = [data[c][i] for c in col[2:]]
                        result = func(*args)
                        self.nuc_del = helper.nuc_del.upper()

                        if result == "empty":
                            break

                        if len(result) != 1:
                            error_reff = error_reff + 1
                            raise ValueError(
                                "HELPER_FUNCTIONS used with vcf_from_snp.py are expected to return a tuple of len 1. Got instead:"
                                + str(result)
                            )
                        var.set_column(vcf_col, result)

                    elif col == "":
                        var.set_column(vcf_col, ".")

                    else:
                        var.set_column(vcf_col, data[col][i])

                var.ref = var.ref.upper()
                var.alt = self.manage_alt(i, var.ref)
                for sample in self.sample_list:
                    var.samples[sample] = self._define_gt(i, var.ref, sample)
                log.debug(f"GT for each sample in variant: {var.samples}")

                if var.ref == "":
                    continue

                if self.del_variant == True:
                    var.ref = self.nuc_del + var.ref

                line = [
                    "chr" + var.chrom,
                    var.pos,
                    var.id,
                    var.ref,
                    var.alt,
                    var.qual,
                ]

                line.append("PASS")
                line.append(".")
                line.append("GT")
                for sample in self.sample_list:
                    line.append(var.samples[sample]["GT"])

                if var.ref == "":
                    raise ValueError("There are no ref for variant", var.chrom, " ", var.pos)

                vcf.write("\t".join(line) + "\n")

        if helper.error_value > 0:
            log.warning(f"{helper.error_value} variant positions not found in ref genome")
