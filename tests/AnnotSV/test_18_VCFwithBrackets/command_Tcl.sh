#!/bin/bash -x

set -eo pipefail

cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"


# The square-bracketed nomenclature ("[[") is used for DEL, DUP, INV, TRA.
# Ex: 
# => Deletion:
# 	1       11099650        .       N       N[1:11101794[     END=11101794;CHR2=1;SVTYPE=DEL;SVLEN=2144
# => Translocation:
# 	1       23224424        .       N       N]6:91004428]     END=91004428;CHR2=6;SVTYPE=TRA;

# AnnotSV must be able to annotate DEL, DUP and INV in their entirety; and the TRA only at the breakend level.



## Checking of the "END" value:
###############################

# -> Check of the end for "TRA" (We should have "end = start" )
# -> Check of the end for the other SV type (end defined with END from INFO)

rm -f "./output/triplesvmerge.CI0_tcl.annotated.tsv" "./output/tutu.CI0"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/triplesvmerge.vcf" -SVinputInfo 1 -outputFile "./output/triplesvmerge.CI0_tcl.annotated.tsv" -includeCI 0 -genomeBuild GRCh37

$cut "./output/triplesvmerge.CI0_tcl.annotated.tsv" "AnnotSV_ID;SV_type;Annotation_mode;SV_start;SV_end;INFO" | grep full > "./output/tutu.CI0"

test=`$ANNOTSV/tests/data/scripts/check_end_values.tcl "./output/tutu.CI0"`

if [ $test == "Ok" ]
then
        echo "Ok"
else
        echo "error 1: END coordinates for bracketed SV not good"
        exit 1
fi




## Vérification avec l'option "-includeCI" à 1 :
################################################

# Si l'on laisse l'option "-includeCI" à 1, end != start+1 dans certains cas pour les TRA (logique !) : On aura $end-$start = $ciend-$cipos+1
# Ce VCF est un merge. On a 2 CIPOS : gridss_CIPOS et brass_CIPOS. Le code utilise le 1er trouvé dans la ligne.
rm -f "./output/triplesvmerge.CI1_tcl.annotated.tsv" "./output/tutu.CI1"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/triplesvmerge.vcf" -SVinputInfo 1 -outputFile "./output/triplesvmerge.CI1_tcl.annotated.tsv" -genomeBuild GRCh37

$cut "./output/triplesvmerge.CI1_tcl.annotated.tsv" "AnnotSV_ID;SV_type;Annotation_mode;SV_start;SV_end;INFO" | grep "TRA" | grep full > "./output/tutu.CI1"

test=`$ANNOTSV/tests/data/scripts/check_includeCI.tcl "./output/tutu.CI1"`

if [ $test == "Ok" ]
then
        echo "Ok"
else
        echo "error 1: END coordinates for bracketed SV not good with -includeCI set to 1"
        exit 1
fi



# GRIDSS output:
# For each of both SV breakends, there are 2 annotation lines:
# 1       11734839        gridss1_4o      C       C[1:11734958[  
# 1       11734958        gridss1_4h      C       ]1:11734839]C  
#################################################################

#$ANNOTSV/bin/AnnotSV -SVinputFile ./input/INPUT_gridss_germline_output.sv.vcf  -SVinputInfo 1 -outputFile ./output/OUTPUT_gridss_germline.sv_tcl.annotated.tsv -genomeBuild GRCh37

# number of variants (SV + SNV/indel) in the input VCF:
#n_VCFinput=`grep -v "^#" ./input/INPUT_gridss_germline_output.sv.vcf | wc -l`
# n_VCFinput => 519

# number of variants in the AnnotSV output:
#n_output=`$cut ./output/OUTPUT_gridss_germline.sv_tcl.annotated.tsv "Annotation_mode" | grep -c full`
# n_output = 283

# number of unnanotated SNV/indel:
#n_nonannotated=`grep -c "< SVminSize " ./output/OUTPUT_gridss_germline.sv_tcl.unannotated.tsv`
# n_nonannotated = 242

#test=`expr $n_output + $n_nonannotated`


echo "ok - Finished"


