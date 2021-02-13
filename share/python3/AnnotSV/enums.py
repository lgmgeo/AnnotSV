from enum import Enum, IntFlag, auto


class BaseEnum(Enum):
    @classmethod
    def metavar(cls):
        """ generates nice metavar strings for typer help """
        return f"[{'|'.join(cls.__members__)}]"

    def __str__(self) -> str:
        return self.name


###


class AnnotationMode(BaseEnum, IntFlag):
    full = auto()
    split = auto()
    both = full | split


class ConfigTypes(str, Enum):
    legacy = "legacy"
    yaml = "yaml"
    json = "json"
    toml = "toml"


class GenomeBuild(str, BaseEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"
    mm9 = "mm9"
    mm10 = "mm10"


class MetricFormat(str, BaseEnum):
    us = "us"
    fr = "fr"


class TranscriptSource(str, BaseEnum):
    RefSeq = "RefSeq"
    ENSEMBL = "ENSEMBL"


class YesNo(BaseEnum):
    yes = True
    no = False
