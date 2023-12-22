#!/bin/bash 

set -eo pipefail

touch "running.flag"

SVinputFile="../test_02_HG00096/input/HG00096.wgs.mergedSV.v8.20130502.svs.genotypes.bed"
SNVindelInputFile="./input/HG00096.chr1.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.withoutSV.filtered.vcf.gz"

$ANNOTSV/bin/AnnotSV -SVinputFile "$SVinputFile" -SVinputInfo 1 -snvIndelFiles "$SNVindelInputFile" -outputFile ./output/HG00096.SV_tcl.annotated.tsv -genomeBuild GRCh37 &> tutu.log 

rm -f "running.flag"


