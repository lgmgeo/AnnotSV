#!/bin/bash

############################################################
# Installing AnnotSV human annotations in a local directory
############################################################

# AIM:
######
# To be used with the "-annotationsDir" option

# CONTEXT:
##########
# To work with bioconda/docker/singularity, AnnotSV couldn't contain the annotations in the recipe (that would make the recipe very large, which is a bad practice in bioconda)
# Users need to download the annotation files once and pass the directory to AnnotSV at runtime with the "-annotationsDir" option.

# USAGE:
########
# INSTALL_annotations.sh "Tag of the AnnotSV version to be cloned"


# CLONE:
########

# cd /path/to/install/annotsv/annotations
mkdir AnnotSV_annotations
cd AnnotSV_annotations

if [ $# -lt 1 ]; then
    echo "Cloning AnnotSV (latest version)"
    git clone https://github.com/lgmgeo/AnnotSV.git
else
    TAG=$1
    echo "Cloning AnnotSV (version $TAG)"
    git clone --branch "$TAG" https://github.com/lgmgeo/AnnotSV.git
fi

cd AnnotSV
make PREFIX=. install
make PREFIX=. install-human-annotation
mv share/AnnotSV/Annotations_Exomiser ..
mv share/AnnotSV/Annotations_Human ..
cd ..
rm -r AnnotSV

