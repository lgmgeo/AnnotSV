#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"

# GRCh37, NM_001037668:
#######################
# NCBI RefSeq genes, curated subset (NM_*, NR_*, NP_* or YP_*):
#   - NM_001037668.1 - chr8:7669242-7673238
# UCSC annotations of RefSeq RNAs (NM_* and NR_*):
#   - NM_001037668 - chr8:7669242-7673238
#   - NM_001037668 - chr8:7353368-7366833


#############################################################################################
## Checking of "-distNearestSS"
#############################################################################################

# Each SV has 2 breakpoints: BND left and BND right


# (Strand -) BND left outside the gene + BND right in the last exon of a gene
#############################################################################

# Input SV: chr8:7669129-7669261

# gr "NM_001037668" $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh37/genes.RefSeq.sorted.bed
# 8       7669241 7673238 -       DEFB107A        NM_001037668    7669329 7673150 7669241,7673080,        7669472,7673238,
# => NM_001037668 : strand = -

# The "SVright" (7669261) is located in the last exon of the transcript NM_001037668 (7669241-7669472)
# The BND is 20 bp after the end of the last exon (3') : 7669261 - 7669241 = 20
#         and 211 bp (7669472 - 7669261)  bp before the start of the exon (last exon => the end of the exon does not correspond to a splice site)
# => we can only calculate the distance in relation to the site in 3' (acceptor)
# => distNearestSS = 211 in 3'

rm -f "./output/test1_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test1.vcf" -outputFile "./output/test1_tcl.annotated.tsv" -genomeBuild GRCh37

# $cut "./output/test1_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Gene_name;Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type"
# AnnotSV_ID      Annotation_mode Gene_name       Tx      Location        Location2       Dist_nearest_SS Nearest_SS_type
# 8_7669129_7669261_DEL_1 full    DEFB107A
# 8_7669129_7669261_DEL_1 split   DEFB107A        NM_001037668    exon2-txEnd     3'UTR   211     3'

annotations=`$cut "./output/test1_tcl.annotated.tsv" "Tx;Dist_nearest_SS;Nearest_SS_type" | grep "NM_001037668"`
if [ `echo $annotations | grep -c "NM_001037668 211 3'"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Dist_nearest_SS and Nearest_SS_type values for NM_001037668 are not correct (test1)"
        exit 1
fi

# SV that overlaps a gene, but with breakpoints outside the gene
################################################################

# We must not have a value for distNearestSS

# gr "bbs10" $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh37/genes.RefSeq.sorted.bed
# 12      76738253        76742195        -       BBS10   NM_024685       76739592        76742138        76738253,76741941,      76741567,76742195,
# => strand = -

# Input SV: 12:76738000-76742200
rm -f "./output/test2_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test2.bed" -outputFile "./output/test2_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37

# $cut "./output/test2_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Gene_name;Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type"
# AnnotSV_ID      Annotation_mode Gene_name       Tx      Location        Location2       Dist_nearest_SS Nearest_SS_type
# 12_76738001_76742200_DEL_1      full    BBS10
# 12_76738001_76742200_DEL_1      split   BBS10   NM_024685       txStart-txEnd   5'UTR-3'UTR

annotations=`$cut "./output/test2_tcl.annotated.tsv" "Tx;Dist_nearest_SS;Nearest_SS_type" | grep "NM_024685"`
if [ $annotations == "NM_024685" ]
then
        echo "Ok"
else
        echo "error 1: Dist_nearest_SS and Nearest_SS_type values for NM_024685 should be empty (test2)"
        exit 1
fi



# SV with breakpoints in a gene "strand -" : 1 BND in the 1st exon, 1 BND in the intron
#######################################################################################

# Input SV: 12:76738280-76741600

# gr "bbs10"  $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh37/genes.RefSeq.sorted.bed
# 12      76738253        76742195        -       BBS10   NM_024685       76739592        76742138        76738253,76741941,      76741567,76742195,
rm -f "./output/test3_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test3.bed" -outputFile "./output/test3_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37

# $cut "./output/test3_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Gene_name;Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type"
# AnnotSV_ID      Annotation_mode Gene_name       Tx      Location        Location2       Dist_nearest_SS Nearest_SS_type
# 12_76738281_76741600_DEL_1      full    BBS10
# 12_76738281_76741600_DEL_1      split   BBS10   NM_024685       intron1-exon2   CDS-3'UTR       33      3'

annotations=`$cut "./output/test3_tcl.annotated.tsv" "Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type" | grep NM_024685`
if [ `echo $annotations | grep -c "NM_024685 intron1-exon2 CDS-3'UTR 33 3'"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Dist_nearest_SS and Nearest_SS_type values for NM_024685 are not correct (test3)"
        exit 1
fi



# SV with breakpoints in a gene containing many exons, strand +
###############################################################

# Input SV: 2:170343580-170349500

# gr "bbs5" $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh37/genes.RefSeq.sorted.bed
# 2       170336003       170363165       +       BBS5    NM_152384       170336063       170361092       170336003,170338760,170343578,170344315,170344496,170349383,170350250,170354136,170355995,170359604,170360812,170360990,    170336122,170338843,170343644,170344365,170344624,170349519,170350346,170354199,170356130,170359688,170360836,170363165,
rm -f "./output/test4_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test4.bed" -outputFile "./output/test4_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37

# $cut "./output/test4_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Gene_name;Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type"
# AnnotSV_ID      Annotation_mode Gene_name       Tx      Location        Location2       Dist_nearest_SS Nearest_SS_type
# 2_170343580_170349500_DEL_1     full    BBS5
# 2_170343580_170349500_DEL_1     split   BBS5    NM_152384       exon3-exon6     CDS     2       3'

annotations=`$cut "./output/test4_tcl.annotated.tsv" "Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type" | grep "NM_152384"`
if [ `echo $annotations | grep -c "NM_152384 exon3-exon6 CDS 2 3'"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Dist_nearest_SS and Nearest_SS_type values for NM_152384 are not correct (test4)"
        exit 1
fi


# SV with a breakpoint located exactly on the 1st base of an exon
#################################################################

# Input SV : 2:170343650-170344315

# gr "bbs5" $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh37/genes.RefSeq.sorted.bed
# 2       170336003       170363165       +       BBS5    NM_152384       170336063       170361092       170336003,170338760,170343578,170344315,170344496,170349383,170350250,170354136,170355995,170359604,170360812,170360990,    170336122,170338843,170343644,170344365,170344624,170349519,170350346,170354199,170356130,170359688,170360836,170363165,
rm -f "./output/test5_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test5.bed" -outputFile "./output/test5_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37

# $cut "./output/test5_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Gene_name;Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type"
# AnnotSV_ID      Annotation_mode Gene_name       Tx      Location        Location2       Dist_nearest_SS Nearest_SS_type
# 2_170343650_170344315_DEL_1     full    BBS5
# 2_170343650_170344315_DEL_1     split   BBS5    NM_152384       intron3-intron3 CDS     0       3'
annotations=`$cut "./output/test5_tcl.annotated.tsv" "Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type" | grep "NM_152384"`
if [ `echo $annotations | grep -c "NM_152384 intron3-intron3 CDS 0 3'"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Dist_nearest_SS and Nearest_SS_type values for NM_152384 are not correct (test5)"
        exit 1
fi



# SV with one breakpoint that matches exactly the end of one exon, and the other breakpoint that matches exactly the start of another exon
##########################################################################################################################################

# Input SV: 11:7716419-7716800 DEL

# grep "NM_198185" $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh37/genes.RefSeq.sorted.bed
# 11      7710667 7728024 -       OVCH2   NM_198185       7711185 7727941 7710667,7711154,7712499,7713132,7716288,7716800,7717231,7718011,7718255,7720296,7721842,7722870,7723262,7723703,7725244,7726111,7727853,    7710833,7711244,7712631,7713226,7716419,7716915,7717257,7718136,7718346,7720320,7722032,7723022,7723358,7723876,7725336,7726221,7728024,

rm -f "./output/test6_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test6.bed" -outputFile "./output/test6_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37

# $cut "./output/test6_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;Gene_name;Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type"
# AnnotSV_ID      Annotation_mode Gene_name       Tx      Location        Location2       Dist_nearest_SS Nearest_SS_type
# 11_7716903_7717219_DEL_1        full    OVCH2
# 11_7716903_7717219_DEL_1        split   OVCH2   NM_198185       intron12-intron12       CDS     0       3'
annotations=`$cut "./output/test6_tcl.annotated.tsv" "Tx;Location;Location2;Dist_nearest_SS;Nearest_SS_type" | grep "NM_198185"`
if [ `echo $annotations | grep -c "NM_198185 intron12-intron12 CDS 0 3'"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Dist_nearest_SS and Nearest_SS_type values for NM_152384 are not correct (test6)"
        exit 1
fi







echo "ok - Finished"



