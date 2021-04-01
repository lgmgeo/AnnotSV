from enum import Enum, IntFlag, auto


class BaseEnum(str, Enum):
    @classmethod
    def metavar(cls):
        """ generates nice metavar strings for typer help """
        return f"[{'|'.join(cls.__members__)}]"

    def __str__(self) -> str:
        return self.name

    @classmethod
    def first(cls):
        return list(cls.__members__.values())[0]


###


class AnnotationMode(BaseEnum):
    full = "full"
    split = "split"
    both = "both"


class ConfigTypes(BaseEnum):
    legacy = "legacy"
    yaml = "yaml"
    json = "json"
    toml = "toml"


class GenomeBuild(BaseEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"
    mm9 = "mm9"
    mm10 = "mm10"


class MetricFormat(BaseEnum):
    us = "us"
    fr = "fr"


class TranscriptSource(BaseEnum):
    RefSeq = "RefSeq"
    ENSEMBL = "ENSEMBL"