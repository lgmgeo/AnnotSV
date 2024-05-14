#!/bin/bash -x

set -eo pipefail

cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


#########################################
# Test with Mouse annotations 
#########################################

rm -f "./output/test_tcl.annotated.tsv" "./output/output.log"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -genomeBuild mm39 -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtBEDcol 4 &> "./output/output.log"

#...output columns annotation (November 08 2021 - 13:34):
#        AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;score;Annotation_mode;CytoBand;Gene_name;Gene_count;Tx;Tx_start;Tx_end;Overlapped_tx_length;Overlapped_CDS_length;Overlapped_CDS_percent;Frameshift;Exon_count;Location;Location2;Dist_nearest_SS;Nearest_SS_type;Intersect_start;Intersect_end;RE_gene;GC_content_left;GC_content_right;Repeat_coord_left;Repeat_type_left;Repeat_coord_right;Repeat_type_right

if [ `grep -c "..AnnotSV is done with the analysis" "./output/output.log"` == 1 ] 
then
        echo "Ok"
else
        echo "error 1: Mouse annotation failed"
        exit 1
fi


# Cytoband annotations:
#######################

# $cut "./output/test_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;CytoBand"
#AnnotSV_ID      Annotation_mode CytoBand
#1_3038297_3156351_DEL_1 full    qA1
#1_94190001_95198035_DUP_1       full    qD
#1_94190001_95198035_DUP_1       split   qD
#...

# list = "CytoBand qA1 qD"
for v in `$cut "./output/test_tcl.annotated.tsv" "CytoBand" | sort -u`
do
        if ! exists_in_list "CytoBand qA1 qD" " " "$v"
        then
                echo "Error with Cytoband annotations"
                echo "error 1"
                exit 1
        fi
done

# Gene_name;Gene_count
######################
# $cut ./output/test_tcl.annotated.tsv "AnnotSV_ID;Annotation_mode;Gene_name;Gene_count"
#--> La ligne full indique 31 gènes


# Tx;Tx_start;Tx_end;Overlapped_tx_length;Overlapped_CDS_length;Overlapped_CDS_percent;Frameshift;Exon_count
############################################################################################################
#$cut ./output/test_tcl.annotated.tsv "AnnotSV_ID;Annotation_mode;Tx;Tx_start;Tx_end;Overlapped_tx_length;Overlapped_CDS_length;Overlapped_CDS_percent;Frameshift;Exon_count"
#AnnotSV_ID      Annotation_mode Tx      Tx_start        Tx_end  Overlapped_tx_length    Overlapped_CDS_length   Overlapped_CDS_percent  Frameshift      Exon_count
#1_3038297_3156351_DEL_1 full
#1_94190001_95198035_DUP_1       full
#1_94190001_95198035_DUP_1       split   NR_138578       94804954        94807341        2387    0       0       no      2
#1_94190001_95198035_DUP_1       split   NR_015525       94830119        94847167        17048   0       0       no      2
#1_94190001_95198035_DUP_1       split   NM_016702       95031816        95041998        10182   1245    100     no      11
#...
#--> 10 colonnes


# Location;Location2;Dist_nearest_SS;Nearest_SS_type
####################################################
#$cut ./output/test_tcl.annotated.tsv "AnnotSV_ID;Annotation_mode;Location;Location2;Dist_nearest_SS;Nearest_SS_type"
#AnnotSV_ID      Annotation_mode Location        Location2       Dist_nearest_SS Nearest_SS_type
#1_3038297_3156351_DEL_1 full
#1_94190001_95198035_DUP_1       full
#1_94190001_95198035_DUP_1       split   txStart-txEnd   UTR
#1_94190001_95198035_DUP_1       split   txStart-txEnd   UTR
#1_94190001_95198035_DUP_1       split   txStart-txEnd   5'UTR-3'UTR
#1_94190001_95198035_DUP_1       split   txStart-txEnd   5'UTR-3'UTR
#...
#1_94190001_95198035_DUP_1       split   exon4-txEnd     CDS-3'UTR       126     3'
#...

# list = "5'UTR-3'UTR CDS-3'UTR Location2 UTR"
for v in `$cut "./output/test_tcl.annotated.tsv" "Location2" | sort -u`
do
        if ! exists_in_list "5'UTR-3'UTR CDS-3'UTR Location2 UTR" " " "$v"
        then
                echo "Error with Location2 annotation"
                echo "error 1"
                exit 1
        fi
done


# Intersect_start;Intersect_end
###############################

#$cut ./output/test_tcl.annotated.tsv "AnnotSV_ID;Annotation_mode;Intersect_start;Intersect_end"
#AnnotSV_ID      Annotation_mode Intersect_start Intersect_end
#1_3038297_3156351_DEL_1 full
#1_94190001_95198035_DUP_1       full
#1_94190001_95198035_DUP_1       split   94804954        94807341
#1_94190001_95198035_DUP_1       split   94830119        94847167
#1_94190001_95198035_DUP_1       split   95031816        95041998
#1_94190001_95198035_DUP_1       split   94766705        94799483
#...


# GC content et repeat
#######################

#$cut ./output/test_tcl.annotated.tsv "AnnotSV_ID;Annotation_mode;GC_content_left;GC_content_right;Repeat_coord_left;Repeat_type_left;Repeat_coord_right;Repeat_type_right" | gr full

#1_3038297_3156351_DEL_1 full    0.290   0.335   1:3038379-3038968       Lx4B
#1_94190001_95198035_DUP_1       full    0.400   0.485                   1:95198045-95198078     (TCC)n


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

# mm10/miRTargetLink_RefSeq_mm10.sorted.bed
# 1       20679055        20679076        mTL_miRNA       Adar;Bdnf;Clcn3;Cpeb1;Eif4e3;Fn1;Fstl1;Fzd7;Gja1;Hdac4;Hhip;Hmgb3;Id1;Id2;Id3;Igsf5;Map4k3;Meox2;Mmd;Mup10;Mup11;Mup12;Mup14;Mup17;Mup19;Nfat5;Notch3;Pax3;Pax7;Pdcd10;Pim1;Pola1;Rarb;Rnf170;Sh3bgrl;Smarcb1;Smarcd2;Spry1;Tgoln1;Timp3;Utrn

# test.RE.bed
# 1       20679054        20689076        DEL


# RE_gene
#########

#$ANNOTSV/bin/AnnotSV -SVinputFile ./input/test.RE.bed -genomeBuild mm10 -SVinputInfo 1 -outputFile ./output/test.RE_tcl.annotated.tsv -REselect1 0 -REselect2 0

#$cut ./output/test.RE_tcl.annotated.tsv "AnnotSV_ID;Annotation_mode;RE_gene" | gr full
#1_20679054_20689076__1  full    Acer1 (RE=mTL_miRNA);Adar (RE=mTL_miRNA);Bdnf (RE=mTL_miRNA);Clcn3 (RE=mTL_miRNA);Cpeb1 (RE=mTL_miRNA);Eif4e3 (RE=mTL_miRNA);Fn1 (RE=mTL_miRNA);Foxl2 (RE=mTL_miRNA);Fstl1 (RE=mTL_miRNA);Fzd7 (RE=mTL_miRNA);Gdf3 (RE=mTL_miRNA);Gja1 (RE=mTL_miRNA);Hdac4 (RE=mTL_miRNA);Hhip (RE=mTL_miRNA);Hmgb3 (RE=mTL_miRNA);Id1 (RE=mTL_miRNA);Id2 (RE=mTL_miRNA);Id3 (RE=mTL_miRNA);Igsf5 (RE=mTL_miRNA);Map4k3 (RE=mTL_miRNA);Meox2 (RE=mTL_miRNA);Mir133b (RE=RefSeq_promoter);Mmd (RE=mTL_miRNA);Mup10 (RE=mTL_miRNA);Mup11 (RE=mTL_miRNA);Mup12 (RE=mTL_miRNA);Mup14 (RE=mTL_miRNA);Mup17 (RE=mTL_miRNA);Mup19 (RE=mTL_miRNA);Nfat5 (RE=mTL_miRNA);Notch3 (RE=mTL_miRNA);Pax3 (RE=mTL_miRNA);Pax7 (RE=mTL_miRNA);Pdcd10 (RE=mTL_miRNA);Pim1 (RE=mTL_miRNA);Pitx3 (RE=mTL_miRNA);Pola1 (RE=mTL_miRNA);Ptbp2 (RE=mTL_miRNA);Rarb (RE=mTL_miRNA);Rnf170 (RE=mTL_miRNA);Sh3bgrl (RE=mTL_miRNA);Smarcb1 (RE=mTL_miRNA);Smarcd2 (RE=mTL_miRNA);Spry1 (RE=mTL_miRNA);Tgoln1 (RE=mTL_miRNA);Timp3 (RE=mTL_miRNA);Tnrc6b (RE=mTL_miRNA);Utrn (RE=mTL_miRNA);Zfp26 (RE=mTL_miRNA)


echo "ok - Finished"

