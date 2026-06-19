#!/bin/bash


# -----------------------------------------------------------------------------
# AIM
# -----------------------------------------------------------------------------
# Download Exomiser and AnnotSV annotations to be used with the "-annotationsDir" option


# -----------------------------------------------------------------------------
# CONTEXT
# -----------------------------------------------------------------------------
# To work with bioconda/docker/singularity, AnnotSV couldn't contain the annotations in the recipe (that would make the recipe very large, which is a bad practice in bioconda)
# Users need to download the annotation files once and pass the directory to AnnotSV at runtime with the "-annotationsDir" option.


# -----------------------------------------------------------------------------
# USAGE
# -----------------------------------------------------------------------------

# Option 1:
# INSTALL_annotations.sh <HUMAN_VERSION> <EXOMISER_VERSION>
#
#   HUMAN_VERSION        Version of AnnotSV human annotations
#   EXOMISER_VERSION     Version of Exomiser phenotype annotations
#
# Example:
#   INSTALL_annotations.sh 3.5 2406
#
#
# Option 2:
# INSTALL_annotations.sh
#
#   In this mode, HUMAN_VERSION and EXOMISER_VERSION are automatically
#   extracted from the AnnotSV Makefile.
#
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# CONFIGURATION (edit before running if needed)
# -----------------------------------------------------------------------------

# ANNOTATIONSDIR: AnnotSV human annotations directory
ANNOTATIONSDIR="./AnnotSV_annotations"

# MAKEFILE: used when HUMAN_VERSION and EXOMISER_VERSION are not provided as arguments
if [[ -z "$1" && -z "$2" ]]; then
	# ANNOTSV: AnnotSV installation directory (must be defined: export ANNOTSV=/path/to/AnnotSV)
	if [ -z "$ANNOTSV" ]; then
	  echo "ERROR: ANNOTSV is not defined"
	  exit 1
	fi

	# AnnotSV Makefile
	MAKEFILE="$ANNOTSV/Makefile"
	
	if [ ! -f "$MAKEFILE" ]; then
		echo "ERROR: Makefile not found at $MAKEFILE"
		echo "Please update the MAKEFILE path in this script (INSTALL_annotations.sh) before running it."
		exit 1
	fi
fi


# -------------------------
# Default HUMAN_VERSION
# -------------------------
if [ -z "$1" ]; then
    HUMAN_VERSION=$(awk -F'= *' '/^HUMAN_VERSION/ {print $2}' "$MAKEFILE")
else
    HUMAN_VERSION="$1"
fi

# -------------------------
# Default EXOMISER_VERSION
# -------------------------
if [ -z "$2" ]; then
	EXOMISER_VERSION=$(awk -F'= *' '/^EXOMISER_VERSION/ {print $2}' "$MAKEFILE")
else
	EXOMISER_VERSION="$2"
fi


mkdir -p $ANNOTATIONSDIR
cd $ANNOTATIONSDIR || exit 1

# ----------------------------------------
# Download AnnotSV human annotations files
# ----------------------------------------
echo ""
echo "Download AnnotSV supporting data files:"
echo ""
curl -C - -LO "https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_${HUMAN_VERSION}.tar.gz"
tar -xf "Annotations_Human_${HUMAN_VERSION}.tar.gz" -C ./
rm -f "Annotations_Human_${HUMAN_VERSION}.tar.gz"

# ---------------------------------------
# Download Exomiser supporting data files
# ---------------------------------------
echo ""
echo "Download Exomiser supporting data files:"
echo ""
curl -C - -LO "https://data.monarchinitiative.org/exomiser/data/${EXOMISER_VERSION}_phenotype.zip"
unzip "${EXOMISER_VERSION}_phenotype.zip" -d "Annotations_Exomiser/${EXOMISER_VERSION}/"
rm -f "${EXOMISER_VERSION}_phenotype.zip"

chmod -R 777 ./Annotations_*



echo "Installation completed. Annotation files have been generated in ${ANNOTATIONSDIR}."


