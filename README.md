# AnnotSV: An integrated tool for Structural Variations annotation and ranking 

## QUICK INSTALLATION
The sources can be cloned to any directory:

cd /’somewhere’/

git clone https://github.com/lgmgeo/AnnotSV.git



Then, the user can easily set the install in the /’somewhere’/AnnotSV/ directory:

cd /'somewhere'/AnnotSV

make PREFIX=. install

make PREFIX=. install-human-annotation

make PREFIX=. install-mouse-annotation

setenv ANNOTSV /’somewhere’/AnnotSV





## TEST

cd /’somewhere’/AnnotSV/share/doc/AnnotSV/Example/

$ANNOTSV/bin/AnnotSV -SVinputFile test.bed -outputFile ./test.annotated.tsv -svtBEDcol 4

