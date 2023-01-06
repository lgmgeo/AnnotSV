# variantconvert

variantconvert is an extendable command-line tool for converting between different file formats used to store genetic variant data. Currently, the following conversions are supported : 

- [AnnotSV](https://lbgi.fr/AnnotSV/) > VCF
- [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) > VCF
- [Arriba](https://github.com/suhrig/arriba/) > VCF
- [DECoN](https://github.com/RahmanTeam/DECoN) > VCF
- BED (including [CANOES](https://github.com/bioinfo-chru-strasbourg/STARK-modules/tree/master/services/structuralvariation/canoes)) > VCF
- [VaRank](https://www.lbgi.fr/VaRank/) (single file or entire folder with varankBatch mode) > VCF

The project is still being developed and maintained. If you encounter any bug or wish to see new features added, please [open an issue](https://github.com/SamuelNicaise/variantconvert/issues).

# Installation

1) Setup an environment with Python >= 3.8. You can use the provided Conda .yml or the Python in the Dockerfile
2) Do the following commands:
```
git clone https://github.com/SamuelNicaise/variantconvert.git
cd variantconvert
pip install -e .
```
3) Change the GENOME["path"] variable in configs/*.json to fit your local system. Some converters will require access to a valid reference genome in fasta format. Alternatively, you can create your own config files in another folder.

# Usage
```
variantconvert --help 
```
Or if you did not use the `pip install` command above:
```
python variantconvert/__main__.py --help
```


____ 

# Section for developers

<details> 
  <summary>Click to read documentation</summary>

## Adding new conversion formats

An intended goal of the project is to make it easy to add new formats to the conversion possibilities. 

Each conversion is described by a JSON config file with the following sections: 
- [GENERAL]
	- skip_rows: how many rows to skip before column indexes
	- unique_variant_id: A list of columns that are needed to uniquely identify a variant. Important for input files where a same variant can be on multiple lines. 

- [VCF_COLUMNS] maps input TSV columns to their corresponding VCF fields. 
  - Add or remove INFO fields at will to customize your output
  - When the equivalence is more complex than 1 input column = 1 VCF field ; you can create advanced HELPER_FUNCTION (explained below).
  
- [COLUMNS_DESCRIPTION]
  - Describes the input tsv columns to write the output VCF header. Column types can be inferred but it is usually safer to define them.

## HELPER_FUNCTION

They're defined in variantconvert/helperfunctions.py and called in your converter's config .json file. 

### To call a HELPER_FUNCTION

Use the following syntax in your .json: 
```
<vcf_field>: ["HELPER_FUNCTION", <function_name>, <tsv column 1>, <tsv column 2>...] # where tsv columns are the TSV fields sent as function input 
```

### To define a HELPER_FUNCTION

1. In HelperFunctions.__init__() , add *<function_name>* to the self.dispatcher dictionary
2. Add a new method in HelperFunctions class named as *<function_name>*, taking as parameters *<tsv column 1>, <tsv column 2>*... in the same order. Then you can use the full power of Python to do any data transformation you wish.

## If customizing a config file is not enough

variantconvert relies on Converter classes that are called by a ConverterFactory depending on the --inputFormat and --outputFormat parameters. 

You can create new Converter classes that will apply different transformations than the existing ones in variantconvert/converters/

They should inherit from the AbstractConverter class and be listed in the ConverterFactory class. That will make them automatically accessible from the command line. 
</details>
