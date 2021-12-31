from __future__ import annotations

import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Tuple

import cyvcf2

if TYPE_CHECKING:
    from annotsv.context import Context


def _strip_chr(val: str):
    return re.sub("chr", "", val, re.I)


@dataclass(init=False)
class Variant:
    __slots__ = ["chrom", "pos", "ref", "alt", "end", "svtype", "svlen", "cipos", "ciend"]
    chrom: str
    pos: int
    ref: str
    alt: str
    end: Optional[int]
    svtype: Optional[str]
    svlen: Optional[int]
    cipos: Optional[Tuple[int, int]]
    ciend: Optional[Tuple[int, int]]

    def __init__(self, var: cyvcf2.Variant):
        self.chrom = _strip_chr(var.CHROM)
        # self.pos = var.POS
        self.pos = var.start + 1
        self.ref = _strip_chr(var.REF)
        if len(var.ALT) == 0:
            self.alt = ""
        else:
            self.alt = _strip_chr(",".join(var.ALT))
        self.svtype = var.INFO.get("SVTYPE")
        self.svlen = var.INFO.get("SVLEN")
        # cyvcf2 uses 0-indexed .end and .start
        self.end = var.INFO.get("END")
        self.cipos = var.INFO.get("CIPOS")
        self.ciend = var.INFO.get("CIEND")

    def __str__(self):
        return self.ident

    @property
    def ident(self):
        return f"{self.chrom}_{self.pos}_{self.end if self.end else ''}_{self.svtype if self.svtype else ''}"


class AnnotationValidator(ABC):
    def __init__(
        self,
        app: Context,
        *,
        label: str,
        downloaded: ResolvedFiles,
        formatted: ResolvedFiles,
    ):
        self._app = app
        self.label = label
        self.downloaded = downloaded
        self.formatted = formatted

    def check(self):
        downloaded = self.downloaded.resolve()
        formatted = self.formatted.resolve()

        success = True
        if not self.downloaded.exists and not self.formatted.exists:
            self._app.log.debug(f"No {self.label} annotation")
            success = False
        elif self.formatted.should_trim():
            self.formatted.keep_last(self.label, self._app)
        elif self.downloaded.exists:
            self.update()
        else:
            self._app.log.debug(f"No new {self.label} annotation to format")
        return success

    @abstractmethod
    def update(self):
        ...

    # annoying, but makes typer checker happy
    def downloaded_path(self):
        return self._get_path(self.downloaded)

    def formatted_path(self):
        return self._get_path(self.formatted)

    def _get_path(self, rf: ResolvedFiles):
        rf_path = rf.resolve()
        if not isinstance(rf_path, Path):
            raise ValueError(f"Got unexpected non-Path downloaded_file: {rf_path!r}")
        return rf_path


class Annotator(ABC):
    @abstractmethod
    def annotate(self, *args, **kwargs):
        ...


class BreakpointAnnotator(Annotator):
    @abstractmethod
    def annotate(self, app: Context, breakpoint_chrom: str, breakpoint_pos: int):
        ...


class SVAnnotator(Annotator):
    @abstractmethod
    def annotate(self, app: Context, sv_chrom: str, sv_start: int, sv_end: int):
        ...


class ResolvedFiles:
    __slots__ = ["dir", "pattern", "extra_patterns", "_members", "_resolved"]
    _members: List[ResolvedFiles]

    def __init__(self, dir_path: Path = None, pattern: str = None, *extra_patterns: str):
        if bool(dir_path) ^ bool(pattern):
            raise ValueError(f"You must set both or neither of dir_path, pattern")
        self.dir = dir_path
        self.pattern = pattern
        self.extra_patterns = list(extra_patterns)
        self._members = []
        self._resolved = None

        if self.extra_patterns and not all("*" in p for p in self.patterns):
            raise ValueError(
                f"Extra patterns can only be used with glob patterns. Got: {self.patterns!r}"
            )

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} dir={self.dir} pattern={self.pattern} _members={len(self._members)}>"

    @property
    def exists(self):
        if self._resolved is None:
            self.resolve()
        assert self._resolved is not None
        if self._members:
            return all([m.exists for m in self._members])
        elif isinstance(self._resolved, list):
            return bool(self._resolved)
        else:
            return self._resolved.exists()

    @property
    def patterns(self):
        if self._members:
            return [m.pattern for m in self._members if m.pattern]
        elif self.pattern:
            return [self.pattern] + self.extra_patterns
        else:
            raise ValueError(f"Cannot return patterns from an empty object")

    def resolve(self, skip_cache: bool = False):
        # cache results for possible expensive globs on slow filesystems
        if not self._resolved or skip_cache:
            if self._members:
                self._resolved = [r.resolve() for r in self._members]
            elif not self.dir or not self.pattern:
                raise ValueError(f"Cannot resolve files with no pattern and/or dir")
            else:
                if "*" in self.pattern:
                    self._resolved = []
                    for subp in self.patterns:
                        self._resolved.extend(list(self.dir.glob(subp)))
                else:
                    self._resolved = self.dir / self.pattern
        return self._resolved

    def add(self, rf: ResolvedFiles):
        if self.dir or self.pattern:
            raise SyntaxError(
                f"Cannot add members to object with existing pattern/dir ({self.pattern!r}/{self.dir!r})"
            )
        self._members.append(rf)

    def should_trim(self):
        if self._members:
            return any(m.should_trim() for m in self._members)
        elif isinstance(self._resolved, list) and len(self._resolved) > 1:
            return True
        else:
            return False

    def keep_last(self, label: str, app: Context):
        if self._members:
            for m in self._members:
                m.keep_last(label, app)
        else:
            assert isinstance(self._resolved, list)
            if len(self._resolved) == 1:
                app.log.debug(f"Not trimming list of single Path")
            else:
                app.keep_last_file(label, self._resolved)
            # if self._resolved is None:
            #     self.resolve()
            # assert self._resolved and self.pattern
            # if isinstance(self._resolved, Path):
            #     app.log.warning(f"Attempted to trim {label} on a single Path ({self._resolved})")
            # else:
            #         app.log.info(
            #             f"Found multiple {label} files, only keeping last: {self._resolved}"
            #         )
            #         for fpath in self._resolved[:-1]:
            #             assert isinstance(
            #                 fpath, Path
            #             ), f"Attempted to prune a list of non-Path objects: {self._resolved}"
            #             new_name = f"{fpath.name}.notused"
            #             fpath.rename(new_name)
