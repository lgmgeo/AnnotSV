#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"
SVinputBEDfile="./input/HG00096.wgs.mergedSV.v8.20130502.svs.genotypes.bed"
SVinputVCFfile="./input/HG00096.wgs.mergedSV.v8.20130502.svs.genotypes.vcf"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


# SVinputFile = BED
pipenv run annotsv --SVinputFile $SVinputBEDfile --svtBEDcol 4 --SVinputInfo 1 --outputFile ./output/HG00096_python_1.annotated.tsv --genomeBuild GRCh37
for v in `$cut ./output/HG00096_python_1.annotated.tsv "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with SVinputFile = BED"
                echo "error 1"
                exit 1
        fi
done

# SVinputFile = VCF
pipenv run annotsv --SVinputFile $SVinputVCFfile --SVinputInfo 1 --outputFile ./output/HG00096_python_2.annotated.tsv --genomeBuild GRCh37
for v in `$cut ./output/HG00096_python_2.annotated.tsv "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with SVinputFile = VCF"
                echo "error 1"
                exit 1
        fi
done

# SVinputFile = zipped VCF
gzip $SVinputVCFfile
pipenv run annotsv --SVinputFile $SVinputVCFfile.gz --SVinputInfo 1 --outputFile ./output/HG00096_python_3.annotated.tsv --genomeBuild GRCh37
for v in `$cut ./output/HG00096_python_3.annotated.tsv "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with SVinputFile = zipped VCF"
                echo "error 1"
                exit 1
        fi
done
gunzip $SVinputVCFfile.gz

# -SVinputInfo 0
pipenv run annotsv --SVinputFile $SVinputBEDfile --svtBEDcol 4 --SVinputInfo 0 --outputFile ./output/HG00096_python_4.annotated.tsv --genomeBuild GRCh37
for v in `$cut ./output/HG00096_python_4.annotated.tsv "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with -SVinputInfo 0"
                echo "error 1"
                exit 1
        fi
done

# -metrics fr
pipenv run annotsv --SVinputFile $SVinputBEDfile --svtBEDcol 4 --outputFile ./output/HG00096_python_5.annotated.tsv --metrics fr --genomeBuild GRCh37
for v in `$cut ./output/HG00096_python_5.annotated.tsv "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with -metrics fr"
                echo "error 1"
                exit 1
        fi
done


# -reciprocal 1
pipenv run annotsv --SVinputFile $SVinputBEDfile --svtBEDcol 4 --SVinputInfo 1 --reciprocal 1 --outputFile ./output/HG00096_python_6.annotated.tsv --genomeBuild GRCh37
for v in `$cut ./output/HG00096_python_6.annotated.tsv "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with -reciprocal 1"
                echo "error 1"
                exit 1
        fi
done


echo "ok - Finished"

