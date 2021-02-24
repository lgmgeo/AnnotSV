#!/usr/bin/env python3

import re
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import List, Optional

import typer

from . import constants
from .config import load_config
from .enums import AnnotationMode, GenomeBuild, MetricFormat, TranscriptSource

### validation / helper funcs


def to_bool(val: str) -> bool:
    try:
        # just use Decimal instead of checking int, float, etc.
        return bool(Decimal(val))
    except InvalidOperation:
        if isinstance(val, str) and val.lower() == "no":
            return False
        return bool(val)


### CLI processing

# init config from file
config = load_config()

# CLI controller
typer_cli = typer.Typer(name="pyAnnotSV")


@typer_cli.command()
def annotsv(
    ctx: typer.Context,
    sv_input_file: Path = typer.Option(
        ...,
        "--SVinputFile",
        help="Path of your VCF or BED input file with SV coordinates",
    ),
    annotations_dir: Optional[Path] = typer.Option(
        str(constants.annotation_dir),
        "--annotationsDir",
        help="Path of the annotations directory",
    ),
    annotation_mode: AnnotationMode = typer.Option(
        config.annotation_mode,
        case_sensitive=False,
        metavar=AnnotationMode.metavar(),
        help="Selection of the type of annotation lines produced by AnnotSV",
    ),
    bcftools: Optional[Path] = typer.Option(
        config.bcftools,
        "--bcftools",
        help="Path of the bcftools local installation",
    ),
    bedtools: Optional[Path] = typer.Option(
        config.bedtools,
        "--bedtools",
        help="Path of the bedtools local installation",
    ),
    candidate_genes_file: Optional[Path] = typer.Option(
        None,
        "--candidateGenesFile",
        dir_okay=False,
        exists=False,
        file_okay=True,
        readable=True,
        resolve_path=True,
        help="Path of a file containing the candidate genes of the user (gene names can be space-separated, tabulation-separated, or line-break-separated) (optional)",
    ),
    candidate_genes_filtering: bool = typer.Option(
        config.candidate_genes_filtering,
        "--candidateGenesFiltering",
        is_flag=False,
        callback=to_bool,
        metavar="0|1",
        help="To select only the SV annotations ('split' and 'full') overlapping a gene from the 'candidateGenesFile'",
    ),
    candidate_snv_indel_files: Optional[str] = typer.Option(
        None,
        "--candidateSnvIndelFiles",
        help="Path of the filtered VCF input file(s) with SNV/indel coordinates for compound heterozygotes report (optional). Gzipped VCF files are supported as well as regular expression",
    ),
    candidate_snv_indel_samples: List[str] = typer.Option(
        [],
        "--candidateSnvIndelsamples",
        help="To specifiy the sample names from the VCF files defined from the -filtereVCFfiles option. Default: use all samples from the filtered VCF files",
    ),
    external_gene_files: List[Path] = typer.Option(
        [],
        "--externalGeneFiles",
        help="Path of tab separated values file(s) to integrate external gene annotations into the output file. The first line should be a header including a column entitled 'genes'. Gzipped files are supported",
    ),
    genome_build: GenomeBuild = typer.Option(
        GenomeBuild.GRCh37,
        "--genomeBuild",
        case_sensitive=False,
        help="Genome build used",
    ),
    hpo: List[str] = typer.Option(
        [],
        "--hpo",
        help="HPO terms list describing the phenotype of the individual being investigated. Values: use comma, semicolon or space separated class values (e.g.: 'HP:0001156,HP:0001363,HP:0011304')",
    ),
    include_ci: bool = typer.Option(
        config.include_ci,
        "--includeCI",
        is_flag=False,
        callback=to_bool,
        help="To expand the 'start' and 'end' SV positions with the VCF confidence intervals (CIPOS, CIEND) around the breakpoints",
    ),
    metrics: MetricFormat = typer.Option(
        MetricFormat.us,
        "--metrics",
        case_sensitive=False,
        help="Changing numerical values from frequencies to us or fr metrics (e.g. 0.2 or 0,2)",
    ),
    min_total_number: int = typer.Option(
        500,
        "--minTotalNumber",
        min=100,
        max=1000,
        help="Minimum number of individuals tested to consider a benign SV for the ranking",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--outputDir",
        dir_okay=True,
        file_okay=False,
        writable=True,
        resolve_path=True,
        help="Output path name",
    ),
    output_file: Path = typer.Option(
        ...,
        "--outputFile",
        dir_okay=False,
        file_okay=True,
        writable=True,
        resolve_path=True,
        help="Output path and file name",
    ),
    overlap: int = typer.Option(
        100,
        "--overlap",
        min=0,
        max=100,
        help="Minimum overlap (%) between the user features and the annotated SV to be reported. Range values: [0-100]",
    ),
    overwrite: bool = typer.Option(
        config.overwrite,
        "--overwrite",
        is_flag=False,
        callback=to_bool,
        help="To overwrite existing output results",
    ),
    promoter_size: int = typer.Option(
        500,
        "--promoterSize",
        min=1,
        help="Number of bases upstream from the transcription start site",
    ),
    rank_filtering: str = typer.Option(
        "1-5,NA",
        "--rankFiltering",
        help="To select the SV of an user-defined specific class (from 1 to 5). Values: use comma separated class values, or use a dash to denote a range of values ( e.g.: '3,4,5' or '3-5')",
    ),
    reciprocal: bool = typer.Option(
        config.reciprocal,
        "--reciprocal",
        is_flag=False,
        callback=to_bool,
        help="Use of a reciprocal overlap between SV and user features (only for annotations with features overlapping the SV)",
    ),
    re_report: bool = typer.Option(
        config.re_report,
        "--REreport",
        is_flag=False,
        callback=to_bool,
        help="Create a report to link the annotated SV and the overlapped regulatory elements (coordinates and sources)",
    ),
    samplesid_bed_col: int = typer.Option(
        -1,
        "--samplesidBEDcol",
        min=4,
        help="Number of the column reporting the samples ID for which the SV was called (if the input SV file is a BED). Range values: [4-[, default = -1 (value not given) $ (Samples ID should be comma or space separated)",
    ),
    snv_indel_files: str = typer.Option(
        ...,
        "--snvIndelFiles",
        help="Path of the VCF input file(s) with SNV/indel coordinates used for false positive discovery. Use counts of the homozygous and heterozygous variants. Gzipped VCF files are supported as well as regular expression",
    ),
    snv_indel_pass: bool = typer.Option(
        False,
        "--snvIndelPASS",
        callback=to_bool,
        help="Boolean. To only use variants from VCF input files that passed all filters during the calling (FILTER column value equal to PASS) $ Values: 0 (default) or 1",
        # NOTE: change to  true/false flag instead of requiring param value?
    ),
    snv_indel_samples: Optional[str] = typer.Option(
        ...,
        "--snvIndelSamples",
        help="To specifiy the sample names from the VCF files defined from the -vcfFiles option. Default: use all samples from the VCF files",
    ),
    sv_input_info: bool = typer.Option(
        True,
        "--SVinputInfo",
        callback=to_bool,
        help="To extract the additional SV input fields and insert the data in the output file $ Values: 1 (default) or 0",
    ),
    sv_min_size: int = typer.Option(50, "--SVminSize", min=1, help="SV minimum size (in bp)"),
    svt_bed_col: int = typer.Option(
        -1,
        "--svtBEDcol",
        min=4,
        help="Number of the column describing the SV type (DEL, DUP) if the input SV file is a BED. Range values: [4-[",
        # NOTE: use Optional/None instead of -1?
    ),
    tx: TranscriptSource = typer.Option(
        TranscriptSource.RefSeq,
        "--tx",
        case_sensitive=False,
        help="Origin of the transcripts",
    ),
    tx_file: Optional[Path] = typer.Option(
        None,
        "--txFile",
        dir_okay=False,
        exists=True,
        file_okay=True,
        readable=True,
        resolve_path=True,
        help="Path of a file containing a list of preferred genes transcripts to be used in priority during the annotation. (Preferred genes transcripts names should be tab or space separated)",
    ),
):
    typer.echo(f"in annotsv")
    breakpoint()
    pass


if __name__ == "__main__":
    # print license info
    typer_cli()
