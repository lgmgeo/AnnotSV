#!/bin/bash

##########################################
# Basic manual install in local directory
##########################################

# Annotations download
######################

# To be used with the "-annotationsDir" option

# Installation of human annotation
mkdir ./AnnotSV_annotations
cd ./AnnotSV_annotations
curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_latest.tar.gz
tar -xf Annotations_Human_latest.tar.gz -C .
rm -rf Annotations_Human_latest.tar.gz


# Installation of Exomiser data:
curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/2202_hg19.tar.gz
curl -C - -LO https://data.monarchinitiative.org/exomiser/data/2202_phenotype.zip
mkdir -p ./Annotations_Exomiser/2202
tar -xf 2202_hg19.tar.gz -C ./Annotations_Exomiser/2202
unzip 2202_phenotype.zip -d ./Annotations_Exomiser/2202
rm -rf 2202_phenotype.zip
rm -rf 2202_hg19.tar.gz


