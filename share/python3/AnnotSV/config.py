from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional

from pydantic import BaseModel, conint, constr, PositiveInt

from . import constants
from .enums import *
from .util import to_camel

# NOTE: defaults coming from CLI, AnnotSV-config.tcl, default configfile. What should be priority? CLI > configfile > class defaults
# should class defaults be eliminated and moved into config? then lose programmatic access to default values though

default_config_file = constants.install_dir / "etc" / "AnnotSV" / "configfile"


class Config(BaseModel):
    # has defaults in config file
    annotation_mode: AnnotationMode
    candidate_genes_filtering: bool
    genome_build: GenomeBuild
    hpo: List[constr(regex=r"^HP:\d+$")]
    include_ci: bool
    metrics: MetricFormat
    min_total_number: conint(ge=100, le=1000)
    overlap: conint(ge=0, le=100)
    overwrite: bool
    promoter_size: PositiveInt
    rank_filtering: List[constr(regex=r"^[1-5]|NA$")]
    reciprocal: bool
    snv_indel_pass: bool
    sv_input_info: bool
    sv_min_size: PositiveInt

    # optional
    annotations_dir: Optional[Path]
    bcftools: Optional[Path]
    bedtools: Optional[Path]
    candidate_genes_file: Optional[Path]
    candidate_snv_indel_files: Optional[str]
    candidate_snv_indel_samples: Optional[List[str]]
    external_gene_files: Optional[List[Path]]
    outputDir: Optional[Path]
    outputFile: Optional[Path]
    re_report: Optional[bool]
    samplesid_bed_col: Optional[int]
    snv_indel_files: Optional[List[Path]]
    snv_indel_samples: Optional[List[str]]
    svt_bed_col: int
    tx: TranscriptSource
    tx_file: Optional[Path]

    output_columns: List

    class Config:
        alias_generator = to_camel


def load_config(
    config_file: Path = default_config_file, config_type: ConfigTypes = ConfigTypes.legacy
) -> Config:
    raise NotImplemented
