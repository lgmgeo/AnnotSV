#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"


# Check of the following output columns:
########################################
# - Overlapped_CDS_length
# - Overlapped_tx_length


# Check on the FCGR3A gene (NM_000569)
######################################
# 
# grep NM_000569 $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh37/genes.RefSeq.sorted.bed
# 1       161511548       161519608       -       FCGR3A  NM_000569       161512801       161519526       161511548,161514490,161518210,161518800,161519486,      161512989,161514748,161518468,161518821,161519608,

# 1       161511549       161519608       -       FCGR3A  NM_000569       161512801       161519526       161511549,161514490,161518210,161518800,161519486,      161512989,161514748,161518468,161518821,161519608,
# 
# Truth = input/test_tcl.annotated.tsv.truth.sauv


rm -f "./output/test_tcl.annotated.tsv"

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/control_9deletions_FCGR3A_NM_000569.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37

compare=`$cut "./output/test_tcl.annotated.tsv"          "Tx;Overlapped_CDS_length;Overlapped_tx_length"`
truth=`$cut "./input/test_tcl.annotated.tsv.truth.sauv"  "Tx;Overlapped_CDS_length;Overlapped_tx_length"`

if [ "$compare" != "$truth" ]
then
        echo "error 1: not the expected values"
        exit 1
else
	echo "ok - Finished"
fi



