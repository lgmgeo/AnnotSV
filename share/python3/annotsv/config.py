import re
import shutil
from pathlib import Path
from typing import List, Optional

from pydantic import BaseModel, Extra, PositiveInt, Field, validator

from annotsv.constants import install_dir
from annotsv.enums import *
from annotsv.util import from_camel, to_camel

# NOTE:
#   config priority: CLI > configfile > class defaults

default_config_file = install_dir / "etc" / "AnnotSV" / "configfile"
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


class Config(BaseModel):
    # has defaults in config file
    annotation_mode: AnnotationMode
    candidate_genes_filtering: bool
    genome_build: GenomeBuild
    hpo: Optional[str]
    include_ci: bool
    metrics: MetricFormat
    min_total_number: int = Field(..., min=100, max=1000)
    overlap: int = Field(..., min=0, max=100)
    overwrite: bool
    promoter_size: PositiveInt
    rank_filtering: str = Field(..., regex=f"^{RANK_PATTERN}(?:,{RANK_PATTERN})*$")
    reciprocal: bool
    snv_indel_pass: bool
    sv_input_info: bool
    sv_min_size: PositiveInt

    # optional
    annotations_dir: Path = Path(__file__).parents[1] / "AnnotSV"
    bcftools: Path = Path(shutil.which("bcftools") or "bcftools")
    bedtools: Path = Path(shutil.which("bedtools") or "bedtools")
    candidate_genes_file: Optional[Path] = None
    candidate_snv_indel_files: Optional[str] = None
    candidate_snv_indel_samples: Optional[List[str]] = None
    external_gene_files: Optional[List[Path]] = None
    outputDir: Optional[Path] = None
    outputFile: Optional[Path] = None
    re_report: Optional[bool] = None
    samplesid_bed_col: Optional[int] = None
    snv_indel_files: Optional[List[Path]] = None
    snv_indel_samples: Optional[List[str]] = None
    svt_bed_col: int
    tx: TranscriptSource
    tx_file: Optional[Path] = None

    output_columns: List[str]

    # @validator("rank_filtering")
    # def validate_rank(cls, v: str) -> str:
    #     rank_pattern = re.compile(r"^([1-5](-[1-5])?)|NA$")
    #     assert all([rank_pattern.search(r) for r in v.split(",")]), f"Invalid ranks in {v}"
    #     return v

    @validator("hpo")
    def validate_hpo(cls, v: Optional[str]):
        if v is None or v == "":
            return None
        assert re.search(f"^{HPO_PATTERN}(?:,{HPO_PATTERN})*$", v), f"Invalid HPO patterns in {v}"
        return v

    # @validator("bcftools", "bedtools")
    # def validate_executable(cls, path_str: str):
    #     bin_path = Path(path_str)
    #     if bin_path.exists():
    #         assert os.access(
    #             bin_path, os.F_OK | os.X_OK
    #         ), f"{bin_path} exists, but is not a file or not executable"
    #     elif bin_path.name == path_str:
    #         which_path = shutil.which(path_str)
    #         assert which_path is not None, f"{path_str} not found"
    #         bin_path = Path(which_path)
    #     else:
    #         assert False, ""
    #     return bin_path.resolve()

    class Config:
        # allows loading/dumping with camelCase, while still keeping snake_case
        alias_generator = to_camel
        # fail on fields that don't exist
        extra = Extra.forbid
        # run validation when changing values on existing config
        validate_assignment = True


def load_config(
    config_file: Path = default_config_file,
    config_type: ConfigTypes = ConfigTypes.LEGACY,
):
    if config_type == ConfigTypes.LEGACY:
        config_dict = _load_legacy_config(config_file)
    else:
        raise NotImplementedError
    config = Config.parse_obj(config_dict)
    return config


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
                elif "$ref" in Config.schema()["properties"][option_name]:
                    # force case-insensitivity in config for enum keys
                    option_value = option_value.upper()
                # invalid options will be caught when creating Config object
                option_config[option_name] = option_value.replace('"', "")
            else:
                raise ValueError(f"unexpected line format: {line}")
    option_config["outputColumns"] = column_config
    return option_config
