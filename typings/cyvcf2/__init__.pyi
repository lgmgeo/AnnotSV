from __future__ import annotations

from pathlib import Path
from typing import (
    Any,
    AnyStr,
    Dict,
    IO,
    Iterable,
    Iterator,
    Mapping,
    MutableMapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)

from numpy import ndarray
from typing_extensions import Literal

Text = AnyStr
Primitives = Union[int, float, bool, AnyStr]


class VCF:
    """
    VCF class holds methods to iterate over and query a VCF.

    Parameters
    ----------
    fname: str
        path to file
    gts012: bool
        if True, then gt_types will be 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN. If False, 3, 2 are flipped.
    lazy: bool
        if True, then don't unpack (parse) the underlying record until needed.
    strict_gt: bool
        if True, then any '.' present in a genotype will classify the corresponding element in the gt_types array as UNKNOWN.
    samples: list
        list of samples to extract from full set in file.
    threads: int
        the number of threads to use including this reader.
    """

    def __init__(
        self,
        fname: Union[str, IO, Path],
        mode: str = "r",
        gts012: bool = False,
        lazy: bool = False,
        strict_gt: bool = False,
        samples: Optional[Sequence[str]] = None,
        threads: Optional[int] = None,
    ) -> None:
        ...

    def __iter__(self) -> "VCF":
        ...

    def __next__(self) -> Variant:
        ...

    def add_filter_to_header(self, adict: Mapping[str, str]) -> "VCF":
        """Add a FILTER line to the VCF header.

        Parameters:
            adict: dict
                dict containing keys for ID, Description."""
        ...

    def add_format_to_header(self, adict: Mapping[str, str]) -> "VCF":
        """Add a FORMAT line to the VCF header.

        Parameters:
            adict: dict
                dict containing keys for ID, Number, Type, Description.
        """
        ...

    def add_info_to_header(self, adict: Mapping[str, str]) -> "VCF":
        """Add a INFO line to the VCF header.

        Parameters:
            adict: dict
                dict containing keys for ID, Number, Type, Description.
        """
        ...

    def add_to_header(self, line: str) -> "VCF":
        """Add a new line to the VCF header

        Parameters:
            line (str) - full vcf header line
        """
        ...

    def contains(self, id: str) -> bool:
        "Check if the given ID is in the header"
        ...

    def gen_variants(
        self,
        sites: Optional[Union[str, Sequence[str]]],
        offset: int = 0,
        each: int = 1,
        call_rate: float = 0.8,
    ) -> Iterable[Tuple[int, Primitives]]:
        ...

    def get_header_type(self, key: str, order: Sequence[int] = [1, 2]) -> Dict[str, str]:
        """
        Extract a field from the VCF header by id.

        Parameters
            key (str) - ID to pull from the header.
            order
        Returns
            dictionary containing header information.
        """
        ...

    def header_iter(self) -> Iterator["HREC"]:
        """Iterate over fields in the HEADER"""
        ...

    @property
    def raw_header(self) -> str:
        """string of the raw header from the VCF"""
        ...

    def relatedness(
        self,
        n_variants: int = 35000,
        gap: int = 30000,
        min_af: float = 0.04,
        max_af: float = 0.8,
        linkage_max: float = 0.2,
        min_depth: int = 8,
    ) -> Mapping[str, Primitives]:
        ...

    def set_index(self, index_path="") -> None:
        ...

    def set_samples(self, samples: Optional[Sequence[AnyStr]]):
        """Set the samples to be pulled from the VCF; this must be called before any iteration.

        Parameters:
            samples (list) - list of samples to extract.
        """
        ...

    @property
    def samples(self) -> Sequence[str]:
        "list of samples pulled from the VCF."
        ...

    @property
    def seqnames(self) -> Sequence[str]:
        "list of chromosomes in the VCF"
        ...


class INFO(MutableMapping):
    """
    INFO is created internally by accessing `Variant.INFO`

    It acts like a dictionary where keys are expected to be in the INFO field of the Variant
    and values are typed according to what is specified in the VCF header

    Items can be deleted with del v.INFO[key] and accessed with v.INFO[key] or v.INFO.get(key)
    """

    ...


class Variant:
    """Variant represents a single VCF Record.
    It is created internally by iterating over a VCF.

    INFO        a dictionary-like field that provides access to the VCF INFO field.
    ALT         the list of alternate alleles.
    CHROM       Chromosome of the variant.
    FILTER      the value of FILTER from the VCF field.
                a value of PASS or '.' in the VCF will give None for this function
    FORMAT      VCF FORMAT field for this variant.
    ID          the value of ID from the VCF field.
    QUAL        the float value of QUAL from the VCF field.
    REF         the reference allele.
    aaf         alternate allele frequency across samples in this VCF.
    call_rate   proportion of samples that were not UNKNOWN.
    end         end of the variant. the INFO field is parsed for SVs."""

    INFO: INFO
    ALT: Sequence[str]
    CHROM: str
    FILTER: Optional[str]
    FORMAT: Mapping[str, Any]
    ID: str
    QUAL: float
    REF: str
    aaf: float
    call_rate: float
    end: int
    POS: int

    def format(
        self,
        field: str,
        vtype: Optional[
            Union[
                type, Literal["Integer"], Literal["Float"], Literal["String"], Literal["Character"]
            ]
        ] = None,
    ) -> ndarray:
        """
        format returns a numpy array for the requested field.
        The numpy array shape will match the requested field. E.g. if the fields has number=3, then the shape will be (n_samples, 3).

        Parameters:
            field (str) - FORMAT field to get the values.
        Return type
            numpy array.
        """
        ...

    @property
    def genotypes(self) -> Sequence[Sequence[Union[int, bool]]]:
        """genotypes returns a list for each sample Indicating the allele and phasing.
        e.g. [0, 1, True] corresponds to 0|1 while [1, 2, False] corresponds to 1/2"""
        ...

    @property
    def gt_alt_depths(self) -> Sequence[int]:
        "get the count of alternate reads as a numpy array."
        ...

    @property
    def gt_alt_freqs(self) -> Sequence[float]:
        "get the freq of alternate reads as a numpy array."
        ...

    @property
    def gt_bases(self) -> Sequence[str]:
        "numpy array indicating the alleles in each sample."
        ...

    @property
    def gt_depths(self) -> Sequence[int]:
        """get the read-depth for each sample as a numpy array."""
        ...

    @property
    def gt_phases(self) -> Sequence[bool]:
        """get a boolean indicating whether each sample is phased as a numpy array."""
        ...

    @property
    def gt_phred_ll_het(self) -> Sequence[Union[int, float]]:
        """get the PL of het for each sample as a numpy array."""
        ...

    @property
    def gt_phred_ll_homalt(self) -> Sequence[Union[int, float]]:
        """get the PL of hom_alt for each sample as a numpy array."""
        ...

    @property
    def gt_phred_ll_homref(self) -> Sequence[Union[int, float]]:
        """get the PL of Hom ref for each sample as a numpy array."""
        ...

    @property
    def gt_quals(self) -> Sequence[float]:
        """get the GQ for each sample as a numpy array."""
        ...

    @property
    def gt_ref_depths(self) -> Sequence[int]:
        """get the count of reference reads as a numpy array."""
        ...

    @property
    def gt_types(self) -> Sequence[int]:
        """gt_types returns a numpy array indicating the type of each sample.
        HOM_REF=0, HET=1. For gts012=True HOM_ALT=2, UNKNOWN=3"""
        ...

    @property
    def is_deletion(self) -> Sequence[bool]:
        """boolean indicating if the variant is a deletion."""
        ...

    @property
    def is_indel(self) -> Sequence[bool]:
        """boolean indicating if the variant is an indel."""
        ...

    @property
    def is_snp(self) -> Sequence[bool]:
        """boolean indicating if the variant is a SNP."""
        ...

    @property
    def is_sv(self) -> Sequence[bool]:
        """boolean indicating if the variant is an SV."""
        ...

    @property
    def is_transition(self) -> Sequence[bool]:
        """boolean indicating if the variant is a transition."""
        ...

    @property
    def num_called(self) -> int:
        """number of samples that were not UNKNOWN."""
        ...

    @property
    def num_het(self) -> int:
        """number heterozygous samples at this variant."""
        ...

    @property
    def num_hom_alt(self) -> int:
        """number homozygous alternate samples at this variant."""
        ...

    @property
    def num_hom_ref(self) -> int:
        """number homozygous reference samples at this variant."""
        ...

    @property
    def num_unknown(self) -> int:
        """number unknown samples at this variant."""
        ...

    @property
    def ploidy(self) -> int:
        """get the ploidy of each sample for the given record."""
        ...

    def set_format(self, name: str, data: ndarray) -> None:
        """set the format field given by name. data must be a numpy array of type float,
        int or string (fixed length ASCII np.bytes_)"""
        ...

    def set_pos(self, pos0: int) -> None:
        "set the POS to the given 0-based position"
        ...

    @property
    def start(self) -> int:
        "0-based start of the variant."
        ...

    @property
    def var_type(self):
        "type of variant (snp/indel/sv)"
        ...


class Writer:
    """
    Writer class makes a VCF Writer.

    Parameters
    ----------
    fname: str
        path to file
    tmpl: VCF
        a template to use to create the output header.
    mode: str
        | Mode to use for writing the file. If ``None`` (default) is given, the mode is
          inferred from the filename extension. If stdout (``"-"``) is provided for ``fname``
          and ``mode`` is left at default, uncompressed VCF will be produced.
        | Valid values are:
        |  - ``"wbu"``: uncompressed BCF
        |  - ``"wb"``: compressed BCF
        |  - ``"wz"``: compressed VCF
        |  - ``"w"``: uncompressed VCF
        | Compression level can also be indicated by adding a single integer to one of
          the compressed modes (e.g. ``"wz4"`` for VCF with compressions level 4).

    Note
    ----
    File extensions ``.bcf`` and ``.bcf.gz`` will both return compressed BCF. If you
    want uncompressed BCF you must explicitly provide the appropriate ``mode``.

    Returns
    -------
    VCF object for iterating and querying.
    """

    name: bytes
    header_written: bool

    def __init__(self, fname: str, tmpl: VCF, mode: str) -> None:
        ...

    def close(self) -> None:
        ...

    @classmethod
    def from_string(cls, fname: AnyStr, header_string: AnyStr, mode: str = "w"):
        ...

    def variant_from_string(self, variant_string: str) -> Variant:
        ...

    def write_header(self):
        ...

    def write_record(self, var: Variant):
        "Write the variant to the writer."
        ...


class HREC:
    @property
    def type(self) -> str:
        """Returns one of: FILTER, INFO, FORMAT, CONTIG, STR, GENERIC"""
        ...

    def __getitem__(self, key: AnyStr) -> str:
        ...

    def info(self, extra: bool = False) -> Mapping[str, Any]:
        """
        return a dict with commonly used stuffs
        """
        ...
