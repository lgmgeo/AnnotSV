<img src="https://raw.githubusercontent.com/SamuelNicaise/variantconvert/master/images/variantconvert_large.png" alt="variantconvert logo"/>

The variantconvert module is an extendable command-line tool for converting between different file formats used to store genetic variant data. Currently, the following conversions are supported :

- [AnnotSV](https://lbgi.fr/AnnotSV/) > VCF
- [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) > VCF
- [Arriba](https://github.com/suhrig/arriba/) > VCF
- [DECoN](https://github.com/RahmanTeam/DECoN) > VCF
- BED (including [CANOES](https://github.com/bioinfo-chru-strasbourg/STARK-modules/tree/master/services/structuralvariation/canoes)) > VCF
- [VaRank](https://www.lbgi.fr/VaRank/) (single file or entire folder with varankBatch mode) > VCF
- [BEDPE](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format) > VCF
- Illumina microarray > VCF

The project is still being developed and maintained. If you encounter any bug or wish to see new features added, please [open an issue](https://github.com/SamuelNicaise/variantconvert/issues).

# Installation

1) Setup an environment with Python >= 3.8. You can use the provided Conda .yml or the Python in the Dockerfile
2) Do the following commands:

```bash
git clone https://github.com/SamuelNicaise/variantconvert.git
cd variantconvert
pip install -e .
```

3) Create a config folder

```bash
# Create a config folder in the default location for your OS.  
variantconvert init

# To see where the folder is: 
variantconvert config --show_config_dir
# All variantconvert commands automatically find this folder without writing the full path (see example in step 4). 

# If you wish, you can move config files to another location and use it instead. 
```

4) Change the GENOME["path"] variable in configs/*.json to fit your local system.

```bash
# Do this for each genome you need to use
# GRCh37
variantconvert config -c GRCh37/* --set GENOME.path=/path/to/your/local/GRCh37.fa --fill_genome_header
# GRCh38
variantconvert config -c GRCh38/* --set GENOME.path=/path/to/your/local/GRCh38.fa --fill_genome_header
# hg19
variantconvert config -c hg19/* --set GENOME.path=/path/to/your/local/hg19.fa --fill_genome_header
# hs1 (T2T)
variantconvert config -c hs1/* --set GENOME.path=/path/to/your/local/hs1.fa --fill_genome_header
# You're ready to use variantconvert !
```

```bash
# You can create your own configs, for example to use other genomes
# Let's create a folder for mm10 (Mus musculus)
cp -r path/to/config/hg19 path/to/config/mm10
variantconvert config -c path/to/config/mm10/* --set GENOME.assembly=mm10 GENOME.path=/path/to/mm10.fa --fill_genome_header
```

Indeed, some converters require a reference genome in fasta format. This is to fill in the VCF "REF" column in cases where we only have the position without the reference base. This implies when converting to a VCF file, you should always use the genome on which the variant was called.

You can create your own config files to customize not only the genome, but also output columns (see the Developers section).

If you installed from source and configured the hg19 genome, you can test that variantconvert is properly installed with the following commands:

```bash
cd <this_repository>
pip install -e .[dev]
pytest
```

# Usage

```bash
variantconvert --help 
```

Or if you did not use the `pip install` command above:

```bash
python src/variantconvert/__main__.py --help
```

Example of a common use case: convert a STAR-Fusion output file to a VCF.

```bash
variantconvert convert -i star-fusions.tsv -o output.vcf -c hg19/starfusion.json
```

<center>

| Conversion  | Default config |
|---|---|
| STAR-Fusion > VCF | starfusion.json  |
| Arriba > VCF | arriba.json  |
| DECoN > VCF | decon.json  |
| BED/CANOES > VCF | canoes_bed.json  |
| BEDPE > VCF | bedpe.json  |
| Illumina microarray > VCF | snp.json  |
|  VaRank > VCF |  varank.json |
|  AnnotSV from BED > VCF | annotsv3_from_bed.json  |
| AnnotSV from VCF > VCF | annotsv3_from_vcf.json  |

</center>

<details>
  <summary>Usage for versions < 2.0.0</summary>

In older versions, input and output format also had to be specified in command line args. Today this is included in config files.

Example of a common use case: convert a STAR-Fusion output file to a VCF.

```txt
variantconvert convert -i star-fusions.tsv -o output.vcf -fi breakpoints -fo vcf -c hg19/starfusion.json
```

<center>

| Conversion  | -fi  (input format) | -fo (output format) | Default config  |
|---|---|---|---|
| STAR-Fusion > VCF  | breakpoints  | vcf  | starfusion.json  |
| Arriba > VCF  | breakpoints  | vcf  | arriba.json  |
| DECoN > VCF  | tsv  | vcf  | decon.json  |
| BED/CANOES > VCF  | tsv  | vcf  | canoes_bed.json  |
| BEDPE > VCF  | bedpe  | vcf  | bedpe.json  |
| Illumina microarray > VCF  | snp  | vcf  | snp.json  |
|  VaRank > VCF | varank  | vcf  |  varank.json |
|  AnnotSV from BED > VCF | annotsv  | vcf  | annotsv3_from_bed.json  |
| AnnotSV from VCF > VCF  | annotsv  | vcf  | annotsv3_from_vcf.json  |

</center>

</details>

<br/>

# Documentation for AnnotSV users

<details>
  <summary>Click to read documentation</summary>

### Creation of a VCF output file format with AnnotSV

To convert the output format from tsv to VCF, AnnotSV relies on the variantconvert tool.s

The variantconvert module distributed with AnnotSV can be used by setting the `-vcf` option to 1 in the AnnotSV command line.

### Requirements in the AnnotSV command line

Different AnnotSV options are required to access to a VCF output:

- From a "BED" or a "VCF" SV input file:
  - The user needs to define the `-SVinputInfo` option to 1 (to report in the tsv output file the 'ID', 'QUAL', 'FILTER'... fields).
- From a "BED" SV input file:
  - The user needs to define the `-svtBEDcol` option (to report the SV type)
  - The `-samplesidBEDcol` option is highly recommended to use (else, the sample colum will be named "NA" (Non Attributed))  

### Method

Each SV from an AnnotSV tsv file is represented with 2 types of lines:

- An annotation on the "full" length of the SV. Every SV are reported, even those not covering a gene.
- An annotation of the SV "split" by gene. This type of annotation gives an opportunity to focus on each gene overlapped by the SV. Thus, when a SV spans over several genes, the output will contain as many annotations lines as genes covered.

#### Example of a duplication overlapping 2 genes (1 full line + 2 split lines in the tsv)

|AnnotSV_ID|SV_chrom|SV_start|SV_end|SV_length|Variant_type|Annotation_mode|Gene_name|DDD_HI_percent|
|---|---|---|---|---|---|---|---|---|
|10_46976157_47590995_1|10|46976157|47590995|614838|DUP|full|AGAP9;ANTXRLP1|91.07|
|10_46976157_47590995_1|10|46976157|47590995|614838|DUP|split|AGAP9|88.1|
|10_46976157_47590995_1|10|46976157|47590995|614838|DUP|split|ANTXRLP1| |

In the converted VCF, each SV is represented with only 1 line by default (mode: "combined" in JSON config). All the annotations (full & split) are reported in the INFO field.
For one SV, all values from a same tsv output column are merged as lists separated by ",". Consequently, all "," in annotations are replaced with "|". If all values (full and all split lines) are identical, they are merged as one.

The tsv output columns are represented in the VCF in this way:
```txt
#mode=combined
#CHROM	POS	REF	ALT	INFO
chr10	46976157	G	<DUP>	AnnotSV_ID=10_46976157_47590995_1;SV_start=46976157;END=47590995;SVLEN=614838;Annotation_mode=full,split,split;Gene_name=AGAP9|ANTXRLP1,AGAP9,ANTXRLP1;DDD_HI_percent=91.07,88.1,.
```

Warning: the AnnotSV > VCF converter uses VCF 4.2 specification, so spaces are replaced with an "_" in the output VCF.

#### Alternative modes

If using lists is complex for your downstream analysis, other conversion modes are available. Instead of combining all full and split annotations, they can be each represented on one line (mode: full&split), or only "full" annotation can be kept (mode: full). Conversion mode can be changed in the JSON config.

```bash
#switch to "full" mode (only full annotations, one line per variant)
variantconvert config -c hg19/* --set GENERAL.mode=full
#switch to "full&split" mode (full and split annotations on separated lines). Use quotes as shown below if using bash because of the '&' character.
variantconvert config -c hg19/* --set 'GENERAL.mode=full&split'
#switch back to default "combined" mode (full and split annotations combined into one line per variant)
variantconvert config -c hg19/* --set GENERAL.mode=combined
```

### GT warning

If the GT is not given in input, the GT is set to "1/." (using the variantconvert distributed by AnnotSV) or "0/1" (using the github variantconvert) for each SV in the VCF output file. Indeed, the considered SV has been called on at least one allele, but we don’t know the status of the second allele. In any case, the user can change this default value in the variantconvert config files.
</details>

<br/>

# Documentation for developers

<details>
  <summary>Click to read documentation</summary>

## Adding new conversion formats

An intended goal of the project is to make it easy to add new formats to the conversion possibilities.

Each conversion is described by a JSON config file with the following sections:

- [GENERAL]
  - `input_format` and `output_format`: Determine which converter module will be returned by ConverterFactory
  - `skip_rows`: how many rows to skip before column indexes
  - `unique_variant_id`: A list of columns that are needed to uniquely identify a variant. Important for input files where a same variant can be on multiple lines.

- [VCF_COLUMNS] maps input TSV columns to their corresponding VCF fields.
  - Add or remove INFO fields at will to customize your output
  - When the equivalence is more complex than 1 input column = 1 VCF field ; you can create advanced HELPER_FUNCTION (explained below).

- [COLUMNS_DESCRIPTION]
  - Describes the input tsv columns to write the output VCF header. Column types can be inferred but it is usually safer to define them.

## HELPER_FUNCTION

They're defined in variantconvert/helpers and called in your converter's config .json file.

### To call a HELPER_FUNCTION

Use the following syntax in your .json:

```bash
<vcf_field>: ["HELPER_FUNCTION", <function_name>, <tsv column 1>, <tsv column 2>...] 
# where tsv columns are the TSV fields sent as function input 
```

### To define a HELPER_FUNCTION

1. In HelperFunctions.\_\_init\_\_() , add *<function_name>* to the self.dispatcher dictionary
2. Add a new method in HelperFunctions class named as *<function_name>*, taking as parameters *<tsv column 1>, <tsv column 2>*... in the same order. Then you can use the full power of Python to do any data transformation you wish.

## If customizing a config file is not enough

variantconvert relies on Converter classes that are called by a ConverterFactory depending on the --inputFormat and --outputFormat parameters (in config file if version >= 2.0.0)

You can create new Converter classes that will apply different transformations than the existing ones in variantconvert/converters/

They should inherit from the AbstractConverter class and be listed in the ConverterFactory class. That will make them automatically accessible from the command line.
</details>
