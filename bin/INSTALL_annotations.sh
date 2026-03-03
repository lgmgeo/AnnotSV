#!/bin/bash

############################################################
# Installing AnnotSV human annotations in a local directory
############################################################

# AIM:
######
# Download Exomiser and AnnotSV annotations to be used with the "-annotationsDir" option

# CONTEXT:
##########
# To work with bioconda/docker/singularity, AnnotSV couldn't contain the annotations in the recipe (that would make the recipe very large, which is a bad practice in bioconda)
# Users need to download the annotation files once and pass the directory to AnnotSV at runtime with the "-annotationsDir" option.

# USAGE:
########
# INSTALL_annotations.sh "Version of AnnotSV human annotation" "Version of Exomiser phenotype annotations"

########

mkdir AnnotSV_annotations
cd AnnotSV_annotations

# AnnotSV annotations
echo ""
echo "Download AnnotSV supporting data files:"
echo ""
curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_$1.tar.gz
tar -xf Annotations_Human_$1.tar.gz -C ./
rm -rf Annotations_Human_$1.tar.gz

#Exomiser
echo ""
echo "Download Exomiser supporting data files:"
echo ""
curl -C - -LO https://data.monarchinitiative.org/exomiser/data/$2_phenotype.zip
unzip $2_phenotype.zip -d Annotations_Exomiser/$2/
rm -rf 2406_phenotype.zip

chmod -R 777 ./Annotations_*




