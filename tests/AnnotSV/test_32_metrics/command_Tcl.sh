#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


# Option "-metrics": changing numerical values from frequencies to us or fr metrics (e.g. 0.2 or 0,2)
# us: 0.2
# fr: 0,2
# WARNING: If the "-metrics" option is set to "fr", this should not affect the ranking


## Checking: 
## - The "." are replaced with ","
## - The ranking stay the same



# -metrics "us" (default)
#########################
rm -f "./output/output.us_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input.bed" -SVinputInfo 1 -outputFile "./output/output.us_tcl.annotated.tsv" -svtbedcol 4 -genomeBuild GRCh37

# $cut "./output/output.us_tcl.annotated.tsv" "B_gain_AFmax;B_loss_AFmax;B_ins_AFmax;B_inv_AFmax;GC_content_left;DDD_HI_percent;ExAC_synZ;GnomAD_pLI;ExAC_pLI;AnnotSV_ranking_score;AnnotSV_ranking_criteria;ACMG_class"
# B_gain_AFmax    B_loss_AFmax    B_ins_AFmax     B_inv_AFmax     GC_content_left DDD_HI_percent  ExAC_synZ       GnomAD_pLI      ExAC_pLI        AnnotSV_ranking_score   AnnotSV_ranking_criteria        ACMG_class
# 0.01    0.0927          0.1667  0.340   0.27    -2.65362022636097       1.0000e+00      1.0000e+00      0.9     1A (cf Gene_count, RE_gene, +0.00);2C-1 (DMD, +0.90);3A (1 gene, +0.00);5F (+0.00)      4
# 0.01    0.0927          0.1667          0.27    -2.65362022636097       1.0000e+00      1.0000e+00                      full=4
#         0.0237          0.0845  0.435   4.17    0.0512652381981393      9.9443e-01      9.9912e-01      0.15    1A (cf Gene_count, RE_gene, +0.00);2H (DIAPH2, +0.15);3A (1 gene, +0.00);5F (+0.00)     3
#         0.0237          0.0845          4.17    0.0512652381981393      9.9443e-01      9.9912e-01                      full=3


for v in `$cut "./output/output.us_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "3 4 ACMG_class full=3 full=4" " " "$v"
        then
                echo "Error with the ranking (-metrics us)"
                echo "error 1"
                exit 1
        fi
done

annotations=`$cut "./output/output.us_tcl.annotated.tsv" "B_gain_AFmax;B_loss_AFmax;B_ins_AFmax;B_inv_AFmax;GC_content_left;DDD_HI_percent;ExAC_synZ;GnomAD_pLI;ExAC_pLI;AnnotSV_ranking_score;AnnotSV_ranking_criteria;ACMG_class" | grep -v "full="`

if [ `echo $annotations | tr " " "\n" | grep -c "\\."` == 25 ] 
then
        echo "Ok"
else
        echo "error 1: Error with the ranking (-metrics us)"
        exit 1
fi



# -metrics "fr"
###############
rm -f "./output/output.fr_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input.bed" -SVinputInfo 1 -metrics fr -outputFile "./output/output.fr_tcl.annotated.tsv" -svtbedcol 4 -genomeBuild GRCh37

# $cut "./output/output.fr_tcl.annotated.tsv" "B_gain_AFmax;B_loss_AFmax;B_ins_AFmax;B_inv_AFmax;GC_content_left;DDD_HI_percent;ExAC_synZ;GnomAD_pLI;ExAC_pLI;AnnotSV_ranking_score;AnnotSV_ranking_criteria;ACMG_class"
# B_gain_AFmax    B_loss_AFmax    B_ins_AFmax     B_inv_AFmax     GC_content_left DDD_HI_percent  ExAC_synZ       GnomAD_pLI      ExAC_pLI        AnnotSV_ranking_score   AnnotSV_ranking_criteria        ACMG_class
# 0,01    0,0927          0,1667  0,340   0,27    -2,65362022636097       1,0     1,0     0,9     1A (cf Gene_count, RE_gene, +0.00);2C-1 (DMD, +0.90);3A (1 gene, +0.00);5F (+0.00)      4
# 0,01    0,0927          0,1667          0,27    -2,65362022636097       1,0     1,0                     full=4
#         0,0237          0,0845  0,435   4,17    0,0512652381981393      0,99443 0,99912 0,15    1A (cf Gene_count, RE_gene, +0.00);2H (DIAPH2, +0.15);3A (1 gene, +0.00);5F (+0.00)     3
#         0,0237          0,0845          4,17    0,0512652381981393      0,99443 0,99912                 full=3


for v in `$cut "./output/output.fr_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "3 4 ACMG_class full=3 full=4" " " "$v"
        then
                echo "Error with the ranking (-metrics us)"
                echo "error 1"
                exit 1
        fi
done

annotations=`$cut "./output/output.fr_tcl.annotated.tsv" "B_gain_AFmax;B_loss_AFmax;B_ins_AFmax;B_inv_AFmax;GC_content_left;DDD_HI_percent;ExAC_synZ;GnomAD_pLI;ExAC_pLI;AnnotSV_ranking_score;AnnotSV_ranking_criteria;ACMG_class" | grep -v "full="`

if [ `echo $annotations | tr " " "\n" | grep -c ","` == 25 ]
then
        echo "Ok"
else
        echo "error 1: Error with the ranking (-metrics fr)"
        exit 1
fi



echo "ok - Finished"



