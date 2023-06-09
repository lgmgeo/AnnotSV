
<p align="center">
    <img src="images/AnnotSV_logo.png" width="500"><br>
</p>
<div align="center">
    <h1 style="font-weight: bold">An integrated tool for Structural Variations annotation and ranking</h1>
</div>
<br>

# Webserver overview
AnnotSV is an annotation engine designed for annotating and ranking Structural Variations (SV). <br>
<br />
To run AnnotSV online, a user-friendly web server interface is freely available [here](https://lbgi.fr/AnnotSV/runjob).

<br />


![](images/AnnotSV_overview.jpg)

<br />

## Input

AnnotSV supports as well the [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) (Variant Call Format) or the commonly used [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) (Browser Extensible Data) input format to describe the SV to annotate (SV type and coordinates). It allows the program to be easily integrated into any bioinformatics pipeline dedicated to NGS analysis.

## Annotations
AnnotSV compiles functionally, regulatory and clinically **relevant information** and aims at providing annotations useful to i) **interpret SV potential pathogenicity and ii) filter out SV potential false positives.**
<br />
These **relevant annotations** are detailed in the [README](README.AnnotSV_3.3.6.pdf).

## Ranking
In addition, AnnotSV provides on top of the annotations **a ranking score** to assess SV pathogenicity.
This score is an adaptation of the work provided by the joint consensus recommendation of ACMG and ClinGen (Riggs et al., 2020). We especially took attention in scoring as much as possible recessive SV observed in various dataset (NGS, array based...).<br />
**The SV classification is detailed [here](ranking.md).**

## Output
The output file contains the overlaps of the SV with relevant genomic features where the genes refer to RefSeq or ENSEMBL genes (user defined).

Various output formats are available online to visualize your AnnotSV results:
-  a TSV (tab-separated values) file powered by the AnnotSV Annotation Engine
-  a VCF file powered by [variantconvert](https://github.com/SamuelNicaise/variantconvert)
-  a knot HTML file and a knot XLSM file powered by [knotAnnotSV](https://github.com/mobidic/knotAnnotSV)
-  an HTML CIRCOS PLOT file powered by [vcf2circos](https://github.com/bioinfo-chru-strasbourg/vcf2circos)


## Typical use
A typical AnnotSV use would be to first look at the annotation and ranking of each SV as a whole (i.e. “full”) and then focus on the content of that SV. Indeed, **there are 2 types of lines produced by AnnotSV** (cf the “AnnotSV type” output column):

- An annotation on the **“full” length** of the SV:
Every SV are reported, even those not covering a gene. This type of annotation gives an estimate of the SV event itself.

- An annotation of the SV **“split” by gene**:
This type of annotation gives an opportunity to focus on each gene overlapped by the SV. Thus, when a SV spans over several genes, the output will contain as many annotations lines as covered genes (cf example in FAQ). This latter annotation is extremely powerful to shorten the identification of mutation implicated in a specific gene.

The annotations columns available in the output file are detailed in the [README](README.AnnotSV_3.3.6.pdf).
