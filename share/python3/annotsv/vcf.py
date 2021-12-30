from __future__ import annotations

import gzip
import re
from pathlib import Path

from cyvcf2 import VCF

from annotsv.constants import VALID_CHROM, VCF_HEADER
from annotsv.context import Context
from annotsv.enums import SVTypes
from annotsv.schemas import Variant
from annotsv.util import append_file, normalize_sv_type, is_empty_file, ymd_hms


def open_vcf(vcf_file: Path):
    if vcf_file.suffix == ".gz":
        return gzip.open(vcf_file, "rt")
    else:
        return vcf_file.open("rt")


def vcf_annotation(chrom: str, start: int, end: int, sv_type: str):
    ...


def vcf2bed(app: Context):
    bedfile = app.config.output_dir / f"{ymd_hms()}_AnnotSV_inputSVfile.bed"
    bed_out = bedfile.open("wt")

    if "annotated" in app.config.output_file.name:
        unannotated_output = app.config.output_dir / app.config.output_file.name.replace(
            "annotated", "unannotated"
        )
    else:
        unannotated_output = app.config.output_file.with_suffix(".unannotated.tsv")

    # save input VCF header for final output
    app.vcf_header = ["SV_type"]
    if "Samples_ID" in app.config.output_columns:
        app.vcf_header.append("Samples_ID")

    vcf = VCF(app.config.sv_input_file)
    if app.config.sv_input_info:
        app.vcf_header.extend(VCF_HEADER[2:])
    else:
        app.vcf_header.extend([VCF_HEADER[3], VCF_HEADER[4], VCF_HEADER[8]])
    app.vcf_header.extend(vcf.samples)

    # start parsing VCF
    gt_missing = True
    for row in vcf:
        # Consider only the SV (not the SNV/indel)
        ##########################################
        # Example of SV:
        # - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG" (length > 50bp)
        # - Type2: alt="<INS>", "<DEL>", ...
        # - Type3: complex rearrangements with breakends: alt="G]17:1584563]" or alt="G]chr17:1584563]"
        var = Variant(row)

        if var.alt.startswith("<"):
            # Type2
            if var.end is None or var.alt == "<TRA>":
                # Example of a strange VCF format line (translocation):
                # chr1   63705386   N    .     <TRA>   51      PASS    "PRECISE;CT=5to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;SVTYPE=TRA;CHR2=chrX;END=478444;SVLEN=0"        GT      ./.
                # POS = 63705386 --> sur chrom1
                # END = 478444   --> sur chromX
                ## => annotation only of the first breakpoint
                var.end = row.POS + 1

            if var.svtype is None:
                var.svtype = var.alt[1:-1]
        elif re.search(r"[\]\[][^:]+:\d+", var.alt):
            # Type3
            # Mainly, that are the breakends which are annotated with this kind of line.
            # But it can also be a DUP, DEL...:
            # 12  46665455  .  N  ]12:46677091]N  .  PASS   END=46677091;CHR2=12;SVTYPE=DUP;SVLEN=11636
            if var.svtype:
                if normalize_sv_type(var.svtype) is SVTypes.NONE or var.end is None:
                    var.end = row.POS + 1
            else:
                var.end = row.POS + 1
        elif re.search(r"^[ACGTN\.\*]+$", f"{var.ref}{var.alt}"):
            # Type1
            refbis = re.sub(r"[\*\.]", "", var.ref)
            altbis = re.sub(r"[\*\.]", "", var.alt)
            var_length = len(altbis) - len(refbis)
            if abs(var_length) < app.config.sv_min_size:
                # it is an indel
                svid = app.get_id(var.ident, var.ref, var.alt)
                append_file(
                    unannotated_output,
                    f"{svid}: variantLength ({abs(var_length)}) < SVminSize ({app.config.sv_min_size})",
                )
                continue
            elif var_length > 0:
                # insertion
                if var.end is None:
                    var.end = var.pos + 1
                if var.svtype is None:
                    var.svtype = "INS"
            else:
                # deletion
                if var.end is None:
                    var.end = var.pos - var_length
                if var.svtype is None:
                    var.svtype = "DEL"

            if "." in var.ref or "." in var.alt:
                # The GRIDSS author says that a . followed by bases refers to a single breakend where the reads cannot be uniquely mapped back.
                # e.g.: 2       39564894        gridss28_45b    T       .TTCTCTCATAACAAACCATGACATCCAGTCATTTAATACAATATGTCTGGGGTGGCTGGGCCCCTTTTTT 246.24  LOW_QUAL
                # => AnnotSV can not determine the SV length
                var_length = None
            elif refbis.lower() not in altbis.lower() and altbis.lower() not in refbis.lower():
                # Complex SV: AGT>ATTGCATGGACCTGAGTCCCCAAAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCGGGGGGGGGG
                # => AnnotSV can not determine the SV length
                var_length = None
            else:
                # AnnotSV can not check the ref and alt values
                # => AnnotSV can not determine the SV length
                var_length = None

            if var.svlen is None and var_length is not None:
                var.svlen = var_length
        else:
            svid = app.get_id(var.ident, var.ref, var.alt)
            append_file(unannotated_output, f"{svid}: not a known SV format")
            continue

        if var.end < var.pos:
            var.pos, var.end = var.end, var.pos
        if var.pos == var.end:
            var.end += 1

        if app.config.include_ci:
            # Correction of the "start" / "end" SV positions by using the confidence interval around the boundaries:
            ##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
            ##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
            # e.g. CIPOS=-10,10;CIEND=-10,10;  or  CIPOS=-6,5;CIEND=-6,5;  or  CIPOS=-410,410;CIEND=-410,410
            if var.cipos:
                ciposleft = int(var.cipos[0])
                if ciposleft < 0:
                    var.pos += ciposleft

            if var.ciend:
                ciendright = int(var.ciend[1])
                if ciendright > 0:
                    var.end += ciendright

        svid = app.get_id(var.ident, var.ref, var.alt)
        if var.chrom not in VALID_CHROM:
            append_file(unannotated_output, f'{svid}: chromosome "{var.chrom}" unknown')
            continue

        if var.end is None:
            append_file(unannotated_output, f"{svid}: END of the SV not defined")
            continue

        if var.svtype is None or var.svtype == "CNV":
            var.svtype = var.alt

        if "GT" in row.FORMAT:
            gt_missing = False
            samples_id = []
            for idx, gt in enumerate(row.genotypes):
                if any([x > 0 for x in gt if not isinstance(x, bool)]):
                    samples_id.append(vcf.samples[idx])
        else:
            # If the GT is not given, AnnotSV considers that the SV is present in all the samples
            samples_id = vcf.samples

        if var.svlen is not None:
            app.sv_lens[svid] = var.svlen

        # Text to write
        raw_line = str(row).strip("\r\n").split("\t")
        bed_line = [var.chrom, var.pos, var.end, var.svtype]
        if "Samples_ID" in app.config.output_columns:
            bed_line.append(",".join(samples_id))

        if app.config.sv_input_info:
            bed_line.extend(raw_line[2:])
        else:
            bed_line.extend(raw_line[8:])
        bed_out.write("\t".join(str(x) for x in bed_line) + "\n")
    bed_out.flush()
    bed_out.close()

    if app.config.candidate_snv_indel_files and gt_missing:
        app.log.warning(
            "The SV genotype is not indicated in the FORMAT column under the 'GT' field. "
            "Compound heterozygosity analysis won't be processed!"
        )
        app.config.candidate_snv_indel_files = None
        app.config.candidate_snv_indel_samples = None

    if is_empty_file(bedfile):
        print("############################################################################")
        print("No SV to annotate in the SVinputFile - Exit without error.")
        print("############################################################################")
        exit(0)

    return bedfile
