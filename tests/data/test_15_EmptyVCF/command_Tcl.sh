#!/bin/bash -x

set -eo pipefail


# VCF input: 
# => contains header lines
# => do not contain variant

# => Error message in the output:
#    ############################################################################
#    SVinputFile (test.vcf is empty, no SV to annotate - Exit without error.
#    ############################################################################

rm -f "./output/test_tcl.annotated.tsv" "./output/test_tcl.log"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.vcf" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" &> "./output/test_tcl.log"

if [ `grep -c "SVinputFile (./input/test.vcf is empty, no SV to annotate - Exit without error." "./output/test_tcl.log"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: \"SVinputFile (test.vcf is empty, no SV to annotate - Exit without error.\" error message NOT found"
        exit 1
fi


echo "ok - Finished"

