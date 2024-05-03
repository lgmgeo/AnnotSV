#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"



####################################
# Checking of Exomiser
####################################


# Test with 1 SV:
#################

# HPO terms = HP:0001156,HP:0001363,HP:0011304,HP:0010055
# Gene positively scored with these HPO = FGFR2 (Entrez Gene: 2263)

# GRCh37:
rm -f "./output/GRCh37-1SV_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/GRCh37-1SV.vcf" -outputFile "./output/GRCh37-1SV_tcl.annotated.tsv" -genomeBuild GRCh37 -hpo "HP:0001156,HP:0001363,HP:0011304,HP:0010055"

annotations=`$cut "./output/GRCh37-1SV_tcl.annotated.tsv" "Annotation_mode;Gene_name;Exomiser_gene_pheno_score;Human_pheno_evidence;Mouse_pheno_evidence;Fish_pheno_evidence"`
# Annotation_mode Gene_name       Exomiser_gene_pheno_score       Human_pheno_evidence    Mouse_pheno_evidence    Fish_pheno_evidence
# full    FGFR2;WDR11;WDR11-AS1   1.0000
#split   FGFR2   1.0000  Brachydactyly;Broad hallux;Broad metatarsal;Broad thumb;Craniosynostosis;Jackson-Weiss syndrome Brachydactyly;Broad hallux;Broad thumb;Craniosynostosis;abnormal sternum morphology;premature cranial suture closure
# split   WDR11   0.6748  Broad hallux;Kallmann syndrome;Pes cavus        Brachydactyly;Broad hallux;Broad thumb;Craniosynostosis;decreased bone mineralization;syndactyly
# split   WDR11-AS1       0.0000

if [ `echo $annotations | grep -c "Craniosynostosis"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Exomiser does not work (GRCh37-1SV)"
        exit 1
fi


# GRCh38:
# => Exomiser annote les gènes (build independant), on doit avoir exactement le même résultat
rm -f "./output/GRCh38-1SV_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/GRCh38-1SV.vcf" -outputFile "./output/GRCh38-1SV_tcl.annotated.tsv" -genomeBuild GRCh38 -hpo "HP:0001156,HP:0001363,HP:0011304,HP:0010055"

annotations=`$cut "./output/GRCh38-1SV_tcl.annotated.tsv" "Annotation_mode;Gene_name;Exomiser_gene_pheno_score;Human_pheno_evidence;Mouse_pheno_evidence;Fish_pheno_evidence"`
if [ `echo $annotations | grep -c "Craniosynostosis"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Exomiser does not work (GRCh38-1SV)"
        exit 1
fi


# Test with > 6000 SV:
######################

# GRCh37:
# Timing: 20 min
rm -f "./output/GRCh37-multi_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/GRCh37-multiSV.vcf" -outputFile "./output/GRCh37-multi_tcl.annotated.tsv" -genomeBuild GRCh37 -hpo "HP:0001156,HP:0001363,HP:0011304,HP:0010055"

annotations=`$cut "./output/GRCh37-multi_tcl.annotated.tsv" "Annotation_mode;Gene_name;Exomiser_gene_pheno_score;Human_pheno_evidence;Mouse_pheno_evidence;Fish_pheno_evidence" | grep PRDM16`
if [ `echo $annotations | grep -c "Craniosynostosis"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Exomiser does not work (GRCh37-multiSV)"
        exit 1
fi


# GRCh38:
# => Exomiser annote les gènes (build independant), on doit avoir exactement le même résultat
# Timing: 20 min
rm -f "./output/GRCh38-multi_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/GRCh38-multiSV.vcf" -outputFile "./output/GRCh38-multi_tcl.annotated.tsv" -genomeBuild GRCh38 -hpo "HP:0001156,HP:0001363,HP:0011304,HP:0010055"

annotations=`$cut "./output/GRCh38-multi_tcl.annotated.tsv" "Annotation_mode;Gene_name;Exomiser_gene_pheno_score;Human_pheno_evidence;Mouse_pheno_evidence;Fish_pheno_evidence" | grep PRDM16`
if [ `echo $annotations | grep -c "Craniosynostosis"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Exomiser does not work (GRCh38-multiSV)"
        exit 1
fi




echo "ok - Finished"



