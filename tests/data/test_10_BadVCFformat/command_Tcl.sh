#!/bin/bash -x

# The return status of the following AnnotSV command is 2 (failed).
# It's normal and we want to continue the script execution => not using the -e option of pipefail
set -o pipefail



# #################################
# ## 1 - Check d'un VCF sans header
# #################################
# 
# The following error message should be displayed:
# ERROR:
# testVCFwithoutHeader.vcf: no VCF header line (prefixed with "#CHROM"). Check your VCF.
# Exit with error
# 
# 
# #################################################################
# # 2 - Non official format for translocations in an input VCF file
# #################################################################
# 
# # Official VCF format
# #####################
# 
# For TRA, the official description in a VCF is done with BND
# 
# 
# # Bad format example:
# #####################
# 
# # Calling of somatic SV with "novobreak" (analysing data from WGS of cancer)
# chr8    1314        N       .       <TRA>   60      PASS    "PRECISE;CT=5to3;CIPOS=-10,10;CIEND=-10,10;SOMATIC;SVTYPE=TRA;CHR2=chr11;END=870;SVLEN=0"   GT      ./.
# 
# => Translocation from chr8 (1314) to chr11 (870)
# => Start = 1314
# => End = 870
#         start > end, pas sur le même chrom
# 
# => From AnnotSV v 2.2, this non official format is correctly parsed by AnnotSV



# VCF without header:
#####################

rm -f "./output/testVCFwithoutHeader_tcl.annotated.tsv" "./output/testVCFwithoutHeader_tcl.log"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/testVCFwithoutHeader.vcf" -SVinputInfo 1 -outputFile "./output/testVCFwithoutHeader_tcl.annotated.tsv" -svtBEDcol 6 -genomeBuild GRCh37 &> "./output/testVCFwithoutHeader_tcl.log"


if [ `grep -c "no VCF header line" "./output/testVCFwithoutHeader_tcl.log"` != 1 ]
then
        echo "error1: Please check the error message for VCF without header."
        exit 1
fi
echo "VCF without header is correctly checked by AnnotSV."



# Bad TRA VCF format: 
#####################

set -eo pipefail

# All the SV should be annotated
rm -f "./output/testTRA_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/testTRA.vcf" -SVinputInfo 1 -outputFile "./output/testTRA_tcl.annotated.tsv" -genomeBuild GRCh37

nInput=`grep -v "^#" "./input/testTRA.vcf" | wc -l`
nOutput=`grep -v "split" "./output/testTRA_tcl.annotated.tsv" | grep -c full`

if [ $nInput != $nOutput ]
then
        echo "error1: Not all SV from the input file ($nInput) are reported in the output ($nOutput)."
        exit 1		
fi


# All the output lines should have the same length
length=`while read line
do  
	echo  "$line" | tr "\t" "\n" | wc -l  
done < "./output/testTRA_tcl.annotated.tsv" | sort -u | wc -l`

if [ $length != 1 ]
then
	echo "error1: The output lines are not all the same length ($length)"
        exit 1
fi



echo "ok - Finished"

