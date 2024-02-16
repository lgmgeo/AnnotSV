#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


# Check that no column is shifted
#################################

# Check the number of annotated SV
###################################
# => should be the same in the VCF and in the annotated file

# There are 2 paires of square-bracketed notations in ALT:
# 
# # First SV:
# 1       67452229        166_1   N       [1:67452635[N   353.30  .       SVTYPE=BND;STRANDS=--:33;CIPOS=-2,0;CIEND=-10,9;CIPOS95=0,0;CIEND95=0,0;MATEID=166_2;EVENT=166;SU=33;PE=17;SR=16        GT:SU:PE:SR:GQ:
# SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB 1/1:9:5:4:87:0.00:-11,-20,-82:102:93:8:92:7:51:3:0:41:4:0.071   0/1:24:12:12:200:353.30:-44,-9,-59:105:81:23:80:22:41:11:0:39:11:0.22
# 1       67452635        166_2   N       [1:67452229[N   353.30  .       SVTYPE=BND;STRANDS=--:33;CIPOS=-10,9;CIEND=-2,0;CIPOS95=0,0;CIEND95=0,0;MATEID=166_1;EVENT=166;SECONDARY;SU=33;PE=17;SR=16      GT:SU:P
# E:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB 0/0:9:5:4:87:0.00:-11,-20,-82:102:93:8:92:7:51:3:0:41:4:0.071   0/1:24:12:12:200:353.30:-44,-9,-59:105:81:23:80:22:41:11:0:39:11:0.22
# 
# # Second SV:
# 1       186826187       370_1   N       N]1:186826404]  98.88   .       SVTYPE=BND;STRANDS=++:8;CIPOS=-10,9;CIEND=-10,9;CIPOS95=0,0;CIEND95=0,0;MATEID=370_2;EVENT=370;SU=8;PE=0;SR=8   GT:SU:PE:SR:GQ:SQ:GL:DP
# :RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB 0/0:3:0:3:48:0.00:-5,-10,-39:47:43:3:43:3:43:3:0:0:0:0.065      0/1:5:0:5:71:98.88:-11,-1,-8:18:12:5:12:5:12:5:0:0:0:0.29
# 1       186826404       370_2   N       N]1:186826187]  98.88   .       SVTYPE=BND;STRANDS=++:8;CIPOS=-10,9;CIEND=-10,9;CIPOS95=0,0;CIEND95=0,0;MATEID=370_1;EVENT=370;SECONDARY;SU=8;PE=0;SR=8 GT:SU:PE:SR:GQ:
# SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB 0/0:3:0:3:48:0.00:-5,-10,-39:47:43:3:43:3:43:3:0:0:0:0.065      0/1:5:0:5:71:98.88:-11,-1,-8:18:12:5:12:5:12:5:0:0:0:0.29




rm -f "./output/test1_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test1.vcf" -SVinputInfo 0 -outputFile "./output/test1_tcl.annotated.tsv" -genomeBuild GRCh37


# Check that no column is shifted
for v in `$cut "./output/test1_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "error1: Error with Lumpy VCF input file"
                exit 1
        fi
done


# Check the number of annotated SV
nVCF=`grep -v "^#" "./input/test1.vcf" | wc -l`

nTSV=`$cut ."/output/test1_tcl.annotated.tsv" "Annotation_mode" | grep -c "full"`
nUnnanotatedTSV=`wc -l "./output/test1_tcl.unannotated.tsv" | cut -d" " -f 1`
nTSV=`echo $(($nTSV + $nUnnanotatedTSV))`

if [ "$nVCF" != "$nTSV" ]
then
        echo "error 1: not the expected number of annotated SV"
        exit 1
else
	echo "ok - Finished"
fi

