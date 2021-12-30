from __future__ import annotations

from enum import Enum, auto


class BaseEnum(str, Enum):
    @classmethod
    def metavar(cls):
        """generates nice metavar strings for typer help"""
        return f"[{'|'.join(cls.__members__)}]"

    @staticmethod
    def _generate_next_value_(name: str, start, count, last_values):
        return name

    def __str__(self) -> str:
        return self.name

    @classmethod
    def first(cls):
        return list(cls.__members__.values())[0]


###


class AnnotationMode(BaseEnum):
    FULL = auto()
    SPLIT = auto()
    BOTH = auto()


class ConfigTypes(BaseEnum):
    LEGACY = auto()
    YAML = auto()
    JSON = auto()
    TOML = auto()


class GenomeBuild(BaseEnum):
    GRCh37 = auto()
    GRCh38 = auto()
    mm9 = auto()
    mm10 = auto()


class MetricFormat(BaseEnum):
    US = auto()
    FR = auto()


class Organisms(BaseEnum):
    Human = auto()
    Mouse = auto()


class TranscriptSource(BaseEnum):
    REFSEQ = auto()
    ENSEMBL = auto()


class SVTypes(BaseEnum):
    DEL = auto()
    DUP = auto()
    INV = auto()
    INS = auto()
    NONE = auto()
