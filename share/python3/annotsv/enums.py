from enum import Enum, auto


class AutoValueMixin:
    """ Makes auto() set the member value to its name """

    @staticmethod
    def _generate_next_value_(name: str, start, count, last_values):
        return name


class BaseEnum(str, AutoValueMixin, Enum):
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
    full = auto()
    split = auto()
    both = auto()


class ConfigTypes(BaseEnum):
    legacy = auto()
    yaml = auto()
    json = auto()
    toml = auto()


class GenomeBuild(BaseEnum):
    GRCh37 = auto()
    GRCh38 = auto()
    mm9 = auto()
    mm10 = auto()


class MetricFormat(BaseEnum):
    us = auto()
    fr = auto()


class TranscriptSource(BaseEnum):
    RefSeq = auto()
    ENSEMBL = auto()
