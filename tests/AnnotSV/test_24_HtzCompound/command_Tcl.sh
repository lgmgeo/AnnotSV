#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"



###############################################################
## Understanding the data used for TESTING
###############################################################


# In "./input/test1.SV.vcf", ther is only 1 SV (identified in sample-002 but not in  sample-001)
################################################################################################################

# CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample-001      sample-002
# 2       27313898        606     N       <DEL>   59.61   .       SVTYPE=DEL;SVLEN=-127;END=27314025;STRANDS=+-:6;CIPOS=-6,5;CIEND=-6,5;CIPOS95=0,0;CIEND95=0,0;SU=6;PE=0;SR=6    GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB     0/0:2:0:2:117:0.00:-4,-16,-61:70:66:3:66:3:66:1:1:0:0:0.043     0/1:4:0:4:59:59.61:-15,-9,-43:61:52:8:52:8:52:6:2:0:0:0.13


# To find the htz compound, we use a VCF (./input/test1.snvindel.vcf) with an indel present in the SV for each of the samples
#############################################################################################################################

# CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample-001      sample-002
# 2       27313998        606     A       ATC     59.61   .       CIPOS95=0,0;CIEND95=0,0;SU=6;PE=0;SR=6  GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB  0/1:2:0:2:117:0.00:-4,-16,-61:70:66:3:66:3:66:1:1:0:0:0.043        0/1:4:0:4:59:59.61:-15,-9,-43:61:52:8:52:8:52:6:2:0:0:0.13


# => In output, we must therefore have an htz compound only for sample-002 (since there is no SV in sample-001) on the "split" line



###############################################################
## TESTS
###############################################################

# The "Sample_ID" feature allows to better visualize the samples which contain the SV (in our example, only the sample "sample-002" contains the SV)
# "Samples_ID" present in the "./input/configfile"

rm -f "./output/test1_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test1.SV.vcf" -candidateSnvIndelFiles "./input/test1.snvindel.vcf" -outputFile "./output/test1_tcl.annotated.tsv" -genomeBuild GRCh37

# AnnotSV_ID                      Samples_ID      Annotation_mode compound_htz(sample-001)        compound_htz(sample-002)        ACMG_class
# 2_27313892_27314030_DEL_1       sample-002      full                                                                            3
# 2_27313892_27314030_DEL_1       sample-002      split                                           2_27313998                      full=3

annotations=`$cut ./output/test1_tcl.annotated.tsv "AnnotSV_ID;Samples_ID;Annotation_mode;compound_htz(sample-001);compound_htz(sample-002);ACMG_class" | grep full`

if [ `echo $annotations | grep -c "sample-002"` == 1 ] && [ `echo $annotations | grep -c "sample-001"` == 0 ]
then
        echo "Ok"
else
        echo "error 1: Exomiser does not work (GRCh37-1SV)"
        exit 1
fi



# We check the cases for which this analysis cannot be carried out:
# - Case 1: with a BED SV input file (-SVinputFile)
# - Case 2: with a SV VCF input file without "GT" (-SVinputFile)
# - Case 3: with a SNV/indel VCF file without "GT" (-candidateSnvIndelFiles)

# Case 1: with a BED SV input file (-SVinputFile)
#################################################

rm -f "./output/test1_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test1.SV.bed" -candidateSnvIndelFiles "./input/test1.snvindel.vcf" -outputFile "./output/test1_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37 &> ./output/output.log

if [ `grep -c "WARNING: The compound heterozygosity analysis is not available from a BED SV input file" ./output/output.log` == 1 ]
then
        echo "Ok"
else
        echo "error 1: No error message displayed for the use of -candidateSnvIndelFiles with a BED input file"
        exit 1
fi


annotations=`$cut "./output/test1_tcl.annotated.tsv" "AnnotSV_ID;Samples_ID;Annotation_mode;compound_htz(sample-001);compound_htz(sample-002);ACMG_class"`

# AnnotSV_ID      Annotation_mode ACMG_class
# 2_27313898_27314025_DEL_1       full    3
# 2_27313898_27314025_DEL_1       split   full=3

if [ `echo $annotations | grep -c "compound_htz"` == 0 ] 
then
        echo "Ok"
else
        echo "error 1: compound_htz annotations seem to be present. Should not with a BED input file"
        exit 1
fi



# Case 2: with a SV VCF input file without "GT" (-SVinputFile)
##############################################################

rm -f "./output/test1_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test1.SV.withoutGT.vcf" -candidateSnvIndelFiles "./input/test1.snvindel.vcf" -outputFile "./output/test1_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37 &> ./output/output.log

#...WARNING: The SV genotype is not indicated in the FORMAT column under the “GT” field
if [ `grep -c "WARNING: The SV genotype is not indicated in the FORMAT column under the " ./output/output.log` == 1 ]
then
        echo "Ok"
else
        echo "error 1: No error message displayed for the use of -candidateSnvIndelFiles with a SV VCF input file without GT"
        exit 1
fi

annotations=`$cut "./output/test1_tcl.annotated.tsv" "AnnotSV_ID;Samples_ID;Annotation_mode;compound_htz(sample-001);compound_htz(sample-002);ACMG_class"`

# AnnotSV_ID      Samples_ID      Annotation_mode ACMG_class
# 2_27313892_27314030_DEL_1       sample-001,sample-002   full    3
# 2_27313892_27314030_DEL_1       sample-001,sample-002   split   full=3

if [ `echo $annotations | grep -c "compound_htz"` == 0 ]
then
        echo "Ok"
else
        echo "error 1: compound_htz annotations seem to be present. Should not with a SV VCF input file without GT"
        exit 1
fi



# Case 3: with a SNV/indel VCF file without "GT" (-candidateSnvIndelFiles)
##########################################################################
rm -f "./output/test1_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test1.SV.vcf" -candidateSnvIndelFiles "./input/test1.snvindel.withoutGT.vcf" -outputFile "./output/test1_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37 &> ./output/output.log

#...WARNING: The SNV/indel genotype is not indicated in the FORMAT column under the “GT” field
if [ `grep -c "WARNING: The SNV/indel genotype is not indicated in the FORMAT column under the" ./output/output.log` == 1 ]
then
        echo "Ok"
else
        echo "error 1: No error message displayed for the use of -candidateSnvIndelFiles with a SNV/indel VCF file without GT"
        exit 1
fi

annotations=`$cut "./output/test1_tcl.annotated.tsv" "AnnotSV_ID;Samples_ID;Annotation_mode;compound_htz(sample-001);compound_htz(sample-002);ACMG_class"`

# AnnotSV_ID      Samples_ID      Annotation_mode ACMG_class
# 2_27313892_27314030_DEL_1       sample-002      full    3
# 2_27313892_27314030_DEL_1       sample-002      split   full=3

if [ `echo $annotations | grep -c "compound_htz"` == 0 ]
then
        echo "Ok"
else
        echo "error 1: compound_htz annotations seem to be present. Should not with a SNV/indel VCF file without GT"
        exit 1
fi



###############################################################
## TESTS WITH 3 samples
###############################################################


# ./input/3samples.theSV.vcf
############################
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  father  mother  prob
# chr1    934051  DEL00000005     G       <DEL>   .       PASS    IMPRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.8.1;CHR2=chr1;END=934869;PE=10;MAPQ=60;CT=3to5;CIPOS=-72,72;CIEND=-72,72;RDRATIO=0.544828    GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV 0/0:-51.0882,0,-14.6956:147:PASS:49:15:32:0:5:10:0:0    0/1:-9.60325,0,-21.6918:96:PASS:53:18:34:0:6:2:0:0      1/0:0,-2.37679,-28.3685:24:PASS:53:30:26:1:8:0:0:0

# DEL: 1:934051-934869
# => detected in "mother,prob"
# => This SV overlaps only 1 gene: SAMD11


# grep SAMD11 $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh38/genes.RefSeq.sorted.bed
# 1       923922  944574  +       SAMD11  NM_001385640    924431  944153  923922,925921,930154,931038,935771,939039,939271,941143,942135,942409,942558,943252,943697,943907,      924948,926013,930336,931089,935896,939129,939412,941306,942251,942488,943058,943377,943808,944574,

# => SAMD11 (GRCh38) NM_001385640, 1:923922-944574



# Search for htz SNVs that are in the overlapped gene (1:923922-944574) (not only in the gene-DEL intersection)
###############################################################################################################

# Extract from the SNV/indel file corresponding to the location of the gene (1:923922-944574) overlapped by the SV
# (We notice that the order of the samples in the SNV/indel file is not the same as in the SV VCF!)

#		prob    father  mother
#chr1 925036     0|0     0|1     0|1
#chr1 925141     1|1     1|1     1|1
#chr1 926250     0|0     0|1     0|1
#chr1 926428     1|1     1|1     1|1
#chr1 926713     1|1     1|1     1|1
#chr1 926744     1|1     1|1     1|1
#chr1 927003     1|1     1|1     1|1
#chr1 927486     0|0     0|1     0|1
#chr1 927744     1|1     1|1     1|1
#chr1 928622     1|1     1|0     1|0
#chr1 929346     1|1     1|1     1|1
#chr1 929375     1|1     1|1     1|1
#chr1 929377     1|1     1|1     1|1
#chr1 929558     0|0     0|1     0|1
#chr1 929839     1|1     1|0     1|0
#chr1 930939     1|1     1|1     1|1
#chr1 931131     0|0     0|1     0|1
#chr1 931513     0|0     0|1     0|1
#chr1 932204     0|0     0|1     0|1
#chr1 932613     0|0     0|1     0|1
#chr1 932949     1|1     1|1     1|1
#chr1 932976     1|1     1|0     1|0
#chr1 933024     0|0     0|1     0|1
#chr1 933411     0|0     0|1     0|1
#chr1 933511     1|1     1|1     1|1
#chr1 933548     0|0     0|1     0|1
#chr1 934601     0|0     0|1     0|1
#chr1 934758     1|0     0|1     0|1
#chr1 935523     0|0     0|1     0|1
#chr1 935954     0|0     0|1     0|1
#chr1 938178     0|0     0|1     0|1
#chr1 939779     0|0     0|1     0|1
#chr1 940262     0|0     0|0     0|1 
#chr1 940263     0/0     0/1     0/2
#chr1 940296     0|0     0|1     ./.
#chr1 940390     0|0     0|1     0|1
#chr1 940820     0|0     0|0     0|1 
#chr1 941119     0|0     0|1     0|1
#chr1 942335     0|0     0|1     0|1
#chr1 943526     0|0     0|1     0|1
#chr1 944296     0|0	0|1	0|1
#chr1 944307     0|0	0|1	0|1

# =>  1 snv/indel htz chez prob
#    28 snv/indel htz chez father
#    29 snv/indel htz chez mother


# AnnotSV
#########

# AnnotSV must report the SNV/indel htz (not the hom) in the mother and the prob (not in the father who does not have this SV)
rm -f "./output/3samples_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/3samples.theSV.vcf" -candidateSnvIndelFiles "./input/3samples.theSnvindel.vcf" -candidateSnvIndelSamples "prob father mother" -genomeBuild GRCh38 -outputFile "./output/3samples_tcl.annotated.tsv" -includeCI 0

n=`$cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;compound_htz(father)" | grep 934051 | tr ";|\t" "\n" | grep -c "1_"`
# 1_934051_934869_DEL_1
# 1_934051_934869_DEL_1
# => ok, nothing reported in the father
if [ $n == 2 ]
then
        echo "Ok"
else
        echo "error 1: some SNV/indel are reported in the father who does not have the SV"
        exit 1
fi


n=`$cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;compound_htz(mother)" | grep 934051 | tr ";|\t" "\n" | grep -c "1_"`
# 1_934051_934869_DEL_1
# 1_934051_934869_DEL_1   1_925036;1_926250;1_927486;1_928622;1_929558;1_929839;1_931131;1_931513;1_932204;1_932613;1_932976;1_933024;1_933411;1_933548;1_934601;1_934758;1_935523;1_935954;1_938178;1_939779;1_940262;1_940263;1_940390;1_940820;1_941119;1_942335;1_943526;1_944296;1_944307
# => ok, 29 snv/indel reporté chez la mother
if [ $n == 31 ]
then
        echo "Ok"
else
        echo "error 1: check the number of SNV/indel reported in the mother (who have the SV)"
        exit 1
fi

n=`$cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;compound_htz(prob)" | grep 934051 | tr ";|\t" "\n" | grep -c "1_"`
# 1_934051_934869_DEL_1
# 1_934051_934869_DEL_1   1_934758
# => ok, 1 snvindel reporté chez le prob
if [ $n == 3 ]
then
        echo "Ok"
else
        echo "error 1: check the number of SNV/indel reported in the probe (who have the SV)"
        exit 1
fi



# On refait en ne reportant que les snv/indel de chez "prob father"
###################################################################

# AnnotSV must report the SNV/indel htz (not the hom) in the prob (not in the father who does not have this SV, not in the mother absent from the "-candidateSnvIndelSamples" option)
rm -f "./output/3samples_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/3samples.theSV.vcf" -candidateSnvIndelFiles "./input/3samples.theSnvindel.vcf" -candidateSnvIndelSamples "prob;father" -genomeBuild GRCh38 -outputFile "./output/3samples_tcl.annotated.tsv" -includeCI 0

n=`$cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;compound_htz(father)" | grep 934051 | tr ";|\t" "\n" | grep -c "1_"`
if [ $n == 2 ]
then
        echo "Ok"
else
        echo "error 1: some SNV/indel are reported in the father who does not have the SV"
        exit 1
fi


n=`$cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;compound_htz(mother)" | grep 934051 | tr ";|\t" "\n" | grep -c "1_"`
if [ $n == 2 ]
then
        echo "Ok"
else
        echo "error 1: some SNV/indel are reported in the mother who should not be reported"
        exit 1
fi

n=`$cut "./output/3samples_tcl.annotated.tsv" "AnnotSV_ID;compound_htz(prob)" | grep 934051 | tr ";|\t" "\n" | grep -c "1_"`
if [ $n == 3 ]
then
        echo "Ok"
else
        echo "error 1: check the number of SNV/indel reported in the probe (who have the SV)"
        exit 1
fi



echo "ok - Finished"




