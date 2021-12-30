from __future__ import annotations

from annotsv.context import Context

# Check the cytoband annotations files
# Create : .../Annotations_$g_AnnotSV(organism)/AnyOverlap/CytoBand/$g_AnnotSV(genomeBuild)/cytoBand_$g_AnnotSV(genomeBuild).formatted.sorted.bed
# (After the formatting step, the "cytoBand_$g_AnnotSV(genomeBuild).bed" is deleted)
def check_cytoband_file(app: Context):
    label = "cytoBand"
    cytoband_bedfile = app.config.cytoband_dir / f"{label}_{app.config.genome_build}.bed"
    formatted_file = cytoband_bedfile.with_suffix(".formatted.bed")
    formatted_sorted_file = cytoband_bedfile.with_suffix(".formatted.sorted.bed")
    header_file = cytoband_bedfile.with_suffix(".header.tsv")

    if cytoband_bedfile.exists():
        raise NotImplementedError()
    elif formatted_sorted_file.exists() and header_file.exists():
        app.log.debug(f"Enabling {label} annotation")
        app.cytoband_ann = True
    else:
        app.log.debug(f"No {label} annotation, ignoring")


def cytoband_annotation(app: Context, breakpoint_chrom: str, breakpoint_pos: int):
    ...
