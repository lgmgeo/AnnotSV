#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"


###################################################################
## Bug: out of memory
## In the VCF input file, some ref or alt values are very larges
###################################################################

## With AnnotSV "2.2" and before
##    couldn't compile regular expression pattern: out of memory
##    while executing
##    "regexp -nocase "$refbis" $altbis"
##    (procedure "VCFsToBED" line 116)


## => Then bugfix added

# ./input/input_1-big-ref.vcf
# => ref with 43367 bp

rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_1-big-ref.vcf" -outputFile "./output/test_tcl.annotated.tsv" -genomeBuild GRCh37 &> ./output/output.log

if [ `grep -c "AnnotSV is done with the analysis" "./output/output.log"` == 1 ] 
then
        echo "Ok"
else
        echo "error 1: AnnotSV does not finished"
        exit 1

fi




echo "ok - Finished"

