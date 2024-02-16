#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"


######################################################################
# Check of the "-snvIndelFiles" option (help for deletion filtering)
######################################################################

# => Add of the columns "Count_hom(HG00096);Count_htz(HG00096);Count_htz/allHom(HG00096);Count_htz/total(cohort);Count_total(cohort)"
# => Values reported only for deletions
# => Display of ratios only if Count_total(cohort) > 50



## ./input/SV.bed
#################
# 1       3777548 4234731 DEL
# 1       3777558 4234720 DUP

# BED => SV location = ]3777548-4234731]


## ./input/extract-HG00096.vcf
##############################
# cut -f 1,2,4,5,10 extract-HG00096.vcf

# #CHROM  POS     REF     ALT     HG00096
# 1       3777548 G       A       1|0:8
# 1       3809347 T       C       0|1:6
# 1       3812848 C       A,G     1|2:8
# 1       3817173 TG      T       1|0:2
# 1       3833398 CCGA    C       0|1:4
# 1       3833877 C       CA,CAA,CAAAA    1|0:0
# 1       3833945 A       AT,ATTT 1|1:4
# 1       4006388 C       A,T     0|1:6
# 1       4015685 C       A       0|1:2
# 1       4026068 G       A,T     1|1:4
# 1       4048855 ATT     AT,A    1|0:8
# 1       4072209 G       A,T     0|2:6
# 1       4073064 G       C,T     0|1:4
# 1       4086587 C       G       0|1:3
# 1       4121636 G       A,T     1|0:3
# 1       4122566 T       TTCTCTCC,TTCTCTCT       1|1:4
# 1       4128523 C       A,G     0|1:5
# 1       4196626 C       A,G     2|0:7
# 1       4204560 TCTTA   CCTTA,T 0|1:15
# 1       4234731 G       C       0|1:5

# => 20 lines <=> 34 SNV/indel (multiallelic) 

# -> The first SNV is not located in the DEL: 1:3777548 (exactly the same position than the DEL start) 
# -> 33 SNV/indel located in the DEL and 32 in the DUP: 65 in total

# It is to notice that AnnotSV reports only 1 count if 2 SNV/indel are located at the same position (e.g. 1       3812848 C       A,G     1|2:8)
# -> DEL : 16 SNV/indel htz + 3 hom



## AnnotSV analysis :
#####################


# First example:
rm -f "./output/SV_tcl.annotated.tsv" "./output/output.log"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/SV.bed" -svtBEDcol 4 -outputFile "./output/SV_tcl.annotated.tsv" -snvIndelFiles "./input/extract-HG00096.vcf" -genomeBuild GRCh37 -includeCI 0 &> "./output/output.log"
#...parsing of snvIndelFiles for "HG00096" (October 29 2021 - 14:43)
#        ...split multiallelic sites into multiple rows for extract-HG00096.vcf (October 29 2021 - 14:43)
#        ...parsing of extract-HG00096.vcf (October 29 2021 - 14:43)
#                -> 65 SNV/indel located in all the SV
#                -> 65 SNV/indel loaded

if [ `grep -c " 65 SNV/indel loaded" "./output/output.log"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: SNV/indel not correctly loaded"
        exit 1
fi


# $cut "./output/SV_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Count_hom(HG00096);Count_htz(HG00096);Count_htz/allHom(HG00096);Count_htz/total(cohort);Count_total(cohort)"
#
# AnnotSV_ID      Annotation_mode Count_hom(HG00096)      Count_htz(HG00096)      Count_htz/allHom(HG00096)       Count_htz/total(cohort) Count_total(cohort)
# 1_3777548_4234731_DEL_1 full    3       16      NA      NA      19
# 1_3777548_4234731_DEL_1 split   0       2       NA      NA      2
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA      0
# 1_3777548_4234731_DEL_1 split   0       1       NA      NA      1
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA      0
# 1_3777548_4234731_DEL_1 split   0       1       NA      NA      1
# 1_3777558_4234720_DUP_1 full    NA      NA      NA      NA      NA
# 1_3777558_4234720_DUP_1 split   NA      NA      NA      NA      NA
# 1_3777558_4234720_DUP_1 split   NA      NA      NA      NA      NA
# 1_3777558_4234720_DUP_1 split   NA      NA      NA      NA      NA
# 1_3777558_4234720_DUP_1 split   NA      NA      NA      NA      NA
# 1_3777558_4234720_DUP_1 split   NA      NA      NA      NA      NA

if [ `$cut "./output/SV_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Count_hom(HG00096);Count_htz(HG00096);Count_htz/allHom(HG00096);Count_htz/total(cohort);Count_total(cohort)" | grep DEL | grep full | cut -f 3` == 3 ]
then
        echo "Ok"
else
        echo "error 1: SNV/indel not correctly loaded (Count_hom)"
        exit 1
fi

if [ `$cut "./output/SV_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Count_hom(HG00096);Count_htz(HG00096);Count_htz/allHom(HG00096);Count_htz/total(cohort);Count_total(cohort)" | grep DEL | grep full | cut -f 4` == 16 ]
then
        echo "Ok"
else
        echo "error 1: SNV/indel not correctly loaded (Count_htz)"
        exit 1
fi



# Check avec 3 samples (trio)
#############################

rm -f "./output/3samples_tcl.annotated.tsv" "./output/output.log"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/3samples.SV.vcf" -outputFile "./output/3samples_tcl.annotated.tsv" -snvIndelFiles "./input/3samples.extract-HG00096.snvindel.vcf" -genomeBuild GRCh38 -metrics us -snvIndelSamples "prob father mother" &> "./output/output.log"

#...parsing of snvIndelFiles for "father mother prob" (October 29 2021 - 16:20)
#        ...split multiallelic sites into multiple rows for 3samples.extract-HG00096.snvindel.vcf (October 29 2021 - 16:20)
#        ...parsing of 3samples.extract-HG00096.snvindel.vcf (October 29 2021 - 16:20)
#                -> 33 SNV/indel located in all the SV (potentially with redundancy depending on SV overlap)
#                -> 33 SNV/indel loaded

if [ `grep -c " 33 SNV/indel loaded" "./output/output.log"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: SNV/indel not correctly loaded"
        exit 1
fi


# $cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Count_hom(father);Count_htz(father);Count_htz/allHom(father);Count_htz/total(cohort)"
#
# AnnotSV_ID      Annotation_mode Count_hom(father)       Count_htz(father)       Count_htz/allHom(father)        Count_htz/total(cohort)
# 1_3777548_4234731_DEL_1 full    3       16      NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   1       4       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA

if [ `$cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Count_hom(father);Count_htz(father);Count_htz/allHom(father);Count_htz/total(cohort)" | grep DEL | grep full | cut -f 3` == 3 ]
then
        echo "Ok"
else
        echo "error 1: SNV/indel not correctly loaded (Count_hom)"
        exit 1
fi

if [ `$cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Count_hom(father);Count_htz(father);Count_htz/allHom(father);Count_htz/total(cohort)" | grep DEL | grep full | cut -f 4` == 16 ]
then
        echo "Ok"
else
        echo "error 1: SNV/indel not correctly loaded (Count_htz)"
        exit 1
fi


# $cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Count_hom(mother);Count_htz(mother);Count_htz/allHom(mother);Count_htz/total(cohort)"
#
# AnnotSV_ID      Annotation_mode Count_hom(mother)       Count_htz(mother)       Count_htz/allHom(mother)        Count_htz/total(cohort)
# 1_3777548_4234731_DEL_1 full    0       19      NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       5       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA

if [ `$cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Count_hom(mother);Count_htz(mother);Count_htz/allHom(mother);Count_htz/total(cohort)" | grep DEL | grep full | cut -f 4` == 19 ]
then
        echo "Ok"
else
        echo "error 1: SNV/indel not correctly loaded (Count_htz)"
        exit 1
fi


# $cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Count_hom(prob);Count_htz(prob);Count_htz/allHom(prob);Count_htz/total(cohort)"
#
# AnnotSV_ID      Annotation_mode Count_hom(prob) Count_htz(prob) Count_htz/allHom(prob)  Count_htz/total(cohort)
# 1_3777548_4234731_DEL_1 full    0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA
# 1_3777548_4234731_DEL_1 split   0       0       NA      NA




echo "ok - Finished"




