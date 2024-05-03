#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"



# With HPO where Exomiser and Phenogenius quite disagree
########################################################
rm -f "./output/test.GRCh37_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile ./input/test.bed  -outputFile "./output/test.GRCh38_tcl.annotated.tsv" -svtbedcol 5 -genomeBuild GRCh38 -hpo "HP:0001156,HP:0001363,HP:0011304,HP:0010055"

# $cut "./output/test.GRCh38_tcl.annotated.tsv" "Annotation_mode;PhenoGenius_score;PhenoGenius_phenotype;PhenoGenius_specificity;Exomiser_gene_pheno_score;Human_pheno_evidence"
# Annotation_mode PhenoGenius_score       PhenoGenius_phenotype   PhenoGenius_specificity Exomiser_gene_pheno_score       Human_pheno_evidence
# full                    C       0.6955
# split   0.47    Brachydactyly   C       0.6955  1p36 deletion syndrome;Brachydactyly;Broad hallux;Broad thumb;Camptodactyly of finger;Craniosynostosis;Delayed cranial suture closure;Foot polydactyly
# split   -1.0                    0.0000
# ...
annotations=`$cut "./output/test.GRCh38_tcl.annotated.tsv" "Annotation_mode;PhenoGenius_specificity" | grep full`
if [ `echo $annotations | grep -c "C"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: check the PhenoGenius_specificity"
        exit 1
fi


# With other HPO where Exomiser and Phenogenius agree
#####################################################
rm -f "./output/PKD1.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile ./input/PKD1.bed -outputFile "./output/PKD1.annotated.tsv" -svtbedcol 4 -genomeBuild GRCh38 -hpo "HP:0000107,HP:0000108,HP:0001407,HP:0005562"
 
# $cut "./output/PKD1.annotated.tsv" "Annotation_mode;PhenoGenius_score;PhenoGenius_phenotype;PhenoGenius_specificity;Exomiser_gene_pheno_score;Human_pheno_evidence"
# Annotation_mode PhenoGenius_score       PhenoGenius_phenotype   PhenoGenius_specificity Exomiser_gene_pheno_score       Human_pheno_evidence
# full                    A       0.9714
# split   6.82    Hepatic cysts, Renal corticomedullary cysts, Multiple renal cysts, Renal cyst   A       0.9714  Hepatic cysts;Multiple renal cysts;Polycystic kidney disease 1;Polycystic kidney dysplasia;Renal corticomedullary cysts;Renal cyst
# ...
annotations=`$cut "./output/PKD1.annotated.tsv" "Annotation_mode;PhenoGenius_specificity" | grep full`
if [ `echo $annotations | grep -c "A"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: check the PhenoGenius_specificity"
        exit 1
fi



echo "ok - Finished"

