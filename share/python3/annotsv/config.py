from __future__ import annotations

import gzip
import os
import re
import shutil
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Extra, Field, PositiveInt, validator

from annotsv import constants
from annotsv.enums import (
    AnnotationMode,
    ConfigTypes,
    GenomeBuild,
    MetricFormat,
    Organisms,
    TranscriptSource,
)
from annotsv.util import from_camel, to_camel

# NOTE:
#   config priority: CLI > configfile > class defaults

default_config_file = constants.install_dir / "etc" / "AnnotSV" / "configfile"
required_output_cols = (
    "AnnotSV_ID",
    "SV_chrom",
    "SV_start",
    "SV_end",
    "SV_length",
    "SV_type",
    "Annotation_mode",
    "Gene_name",
    "Gene_count",
    "Tx",
    "Tx_start",
    "Tx_end",
    "Overlapped_tx_length",
    "Overlapped_CDS_length",
    "Frameshift",
    "Location",
    "Location2",
    "Dist_nearest_SS",
    "Nearest_SS_type",
    "Intersect_start",
    "Intersect_end",
    "RE_gene",
    "B_gain_source",
    "B_gain_coord",
    "B_loss_source",
    "B_loss_coord",
    "HI",
    "TS",
    "GnomAD_pLI",
    "OMIM_morbid",
    "OMIM_morbid_candidate",
    "AnnotSV_ranking_score",
    "AnnotSV_ranking_criteria",
    "ACMG_class",
)

RANK_PATTERN = r"(?:[1-5](?:-[1-5])?)|NA"
HPO_PATTERN = r"HP:\d+$"
ORGANISM_MAP = {
    GenomeBuild.GRCh37: Organisms.Human,
    GenomeBuild.GRCh38: Organisms.Human,
    GenomeBuild.mm9: Organisms.Mouse,
    GenomeBuild.mm10: Organisms.Mouse,
}


class AnnotSVConfig(BaseModel):
    # from command line
    output_file: Path
    output_dir: Path
    sv_input_file: Path
    snv_indel_files: List[Path]
    snv_indel_samples: List[str]

    # has defaults in config file
    annotation_mode: AnnotationMode
    benign_af: float = Field(..., min=0.001, max=0.1)
    candidate_genes_filtering: bool
    genome_build: GenomeBuild
    hpo: Optional[str] = None
    include_ci: bool
    metrics: MetricFormat
    min_total_number: int = Field(..., min=100, max=1000)
    overlap: int = Field(..., min=0, max=100)
    overwrite: bool
    promoter_size: PositiveInt
    rank_filtering: str = Field(..., regex=f"^{RANK_PATTERN}(?:,{RANK_PATTERN})*$")
    reciprocal: bool
    snvindel_pass: bool
    sv_input_info: bool
    sv_min_size: PositiveInt

    # optional
    annotations_dir: Path = Path(__file__).parents[1] / "AnnotSV"
    bcftools: Path = Path("bcftools")
    bedtools: Path = Path("bedtools")
    candidate_genes_file: Optional[Path] = None
    candidate_snv_indel_files: Optional[str] = None
    candidate_snv_indel_samples: Optional[List[str]] = None
    external_gene_files: Optional[List[Path]] = None
    re_report: bool
    samplesid_bed_col: Optional[int] = None
    svt_bed_col: int
    tx: TranscriptSource
    tx_file: Optional[Path] = None

    output_columns: List[str]

    @property
    def organism(self):
        return ORGANISM_MAP[self.genome_build]

    @property
    def annotation_dir(self):
        return constants.annotation_dir / f"Annotations_{self.organism}"

    @property
    def annotation_root(self):
        # consistency / convenience
        return constants.annotation_dir

    @property
    def benign_dir(self):
        return self.benign_root / self.genome_build.name

    @property
    def benign_root(self):
        return self.annotation_dir / f"SVincludedInFt/BenignSV"

    @property
    def blacklist_dir(self):
        return self.breakpoint_dir / f"ENCODEblacklist/{self.genome_build}"

    @property
    def breakpoint_dir(self):
        return self.annotation_dir / "BreakpointsAnnotations"

    @property
    def clingen_dir(self):
        return self.extann_dir / "ClinGen"

    @property
    def cosmic_dir(self):
        return self.annotation_dir / f"FtIncludedInSV/COSMIC/{self.genome_build}"

    @property
    def cytoband_dir(self):
        return self.annotation_dir / f"AnyOverlap/CytoBand/{self.genome_build}"

    @property
    def exac_dir(self):
        return self.extann_dir / "ExAC"

    @property
    def exomiser_dir(self):
        # for consistency
        return constants.exomiser_dir

    @property
    def extann_dir(self):
        return self.annotation_dir / "Gene-based"

    @property
    def gap_dir(self):
        return self.breakpoint_dir / f"Gap/{self.genome_build}"

    @property
    def gcc_dir(self):
        return self.breakpoint_dir / f"GCcontent/{self.genome_build}"

    @property
    def genes_dir(self):
        return self.annotation_dir / f"Genes/{self.genome_build}"

    @property
    def ncbi_dir(self):
        return self.extann_dir / "NCBIgeneID"

    @property
    def omim_dir(self):
        return self.extann_dir / "OMIM"

    @property
    def pathogenic_snv_dir(self):
        return self.pathogenic_snv_root / self.genome_build.name

    @property
    def pathogenic_snv_root(self):
        return self.annotation_dir / "FtIncludedInSV/PathogenicSNVindel"

    @property
    def pathogenic_sv_dir(self):
        return self.pathogenic_sv_root / self.genome_build.name

    @property
    def pathogenic_sv_root(self):
        return self.annotation_dir / "FtIncludedInSV/PathogenicSV"

    @property
    def promoter_dir(self):
        return self.annotation_dir / f"FtIncludedInSV/Promoter/{self.genome_build}"

    @property
    def reg_elements_dir(self):
        return self.annotation_dir / f"FtIncludedInSV/RegulatoryElements/{self.genome_build}"

    @property
    def repeat_dir(self):
        return self.breakpoint_dir / f"Repeat/{self.genome_build}"

    @property
    def segdup_dir(self):
        return self.breakpoint_dir / f"SegDup/{self.genome_build}"

    @property
    def tad_dir(self):
        return self.annotation_dir / f"FtIncludedInSV/TAD/{self.genome_build}"

    @property
    def user_bed_dir(self):
        return self.annotation_dir / f"Users/{self.genome_build}"

    @validator("hpo")
    def validate_hpo(cls, v: Optional[str]):
        if v is None or v == "":
            return None
        if not re.search(f"^{HPO_PATTERN}(?:,{HPO_PATTERN})*$", v):
            raise ValueError(f"Invalid HPO patterns in {v}")
        return v

    @validator("snv_indel_samples", "candidate_snv_indel_samples")
    def validate_samples(cls, v, values, field, **kwargs):
        v_list: List[str] = []
        if v:
            v_list = re.split(r"[;,|]", v)

        files_field = field.name.replace("samples", "files")
        if values.get(files_field):
            # All samples should be present in one of the VCF files
            # If the no samples are defined, include all samples found in the VCFs
            all_samples = set()
            for fpath in values[files_field]:
                if fpath.suffix == ".gz":
                    open_func = gzip.open
                else:
                    open_func = open

                with open_func(fpath, "rt") as fh:
                    for line in fh:
                        if line.startswith("#CHROM"):
                            all_samples.update(line.rstrip().split("\t")[9:])
                            break

            if len(v_list) == 0:
                v_list = list(all_samples)
            else:
                v_list = list(set(v_list) & all_samples)

        return v_list

    @validator("output_file")
    def validate_output_file(cls, v: Path, values, **kwargs):
        if v.suffix != ".tsv":
            v = Path(f"{v}.tsv")
        return v

    @validator("output_dir")
    def validate_output_dir(cls, v, values, **kwargs):
        if v:
            if (
                values["output_file"]
                and values["output_file"].is_absolute()
                and not str(values["output_file"].resolve()).startswith(str(v.resolve()))
            ):
                raise ValueError(
                    f"Output dir {v} does not match output file {values['output_file']}"
                )
        else:
            v = values["output_file"].resolve().parent
        v.mkdir(parents=True, exist_ok=True)
        return v

    @validator("bcftools", "bedtools")
    def validate_executable(cls, v: Path, field, **kwargs):
        if v.exists():
            if not os.access(v, os.F_OK | os.X_OK):
                raise ValueError(f"{v} exists, but is not a file or not executable")
        else:
            which_path = shutil.which(field.name)
            if which_path is None:
                raise ValueError(f"{field.name} not found")
            v = Path(which_path)
        return v.resolve()

    class Config:
        # allows loading/dumping with camelCase, while still keeping snake_case
        alias_generator = to_camel
        # fail on fields that don't exist
        extra = Extra.forbid
        # run validation when changing values on existing config
        validate_assignment = True


def load_config(
    params: Dict[str, Any],
    *,
    config_file: Path = default_config_file,
    config_type: ConfigTypes = ConfigTypes.LEGACY,
):
    if config_type == ConfigTypes.LEGACY:
        config_defaults = _load_legacy_config(config_file)
    else:
        raise NotImplementedError

    # keys present in params will overwrite values from config file
    return AnnotSVConfig.parse_obj({**config_defaults, **params})


def _load_legacy_config(config_file: Path):
    opt_val_pattern = re.compile(r"^-([a-zA-Z]+):\s+(\S+?)\s*$")
    column_pattern = re.compile(r"[ \t\*:]+$")
    with config_file.open("rt") as input:
        option_config = dict()
        column_config = list(required_output_cols)
        column_mode = False
        for line in input:
            if line.startswith("# AnnotSV Output columns:"):
                # switch to output column mode
                column_mode = True
            elif line.lstrip().startswith("#") or line.strip() == "":
                # ignore empty / comment lines
                continue
            elif column_mode:
                # non-empty line in column mode
                col_name = re.sub(column_pattern, "", line.rstrip())
                if col_name not in column_config:
                    column_config.append(col_name)
            elif line.startswith("-"):
                # non-empty line in option mode
                match = opt_val_pattern.search(line.strip())
                if match is None:
                    raise ValueError(
                        f"Invalid config line, does not match /{opt_val_pattern.pattern}/:\n{line}"
                    )
                option_name, option_value = match.group(1, 2)
                if option_name == "rankFiltering":
                    pass
                elif issubclass(AnnotSVConfig.__fields__[from_camel(option_name)].type_, Enum):
                    # force case-insensitivity in config for enum keys
                    option_value = option_value.upper()
                # invalid options will be caught when creating Config object
                option_config[option_name] = option_value.replace('"', "")
            else:
                raise ValueError(f"unexpected line format: {line}")
    option_config["outputColumns"] = column_config
    return option_config
