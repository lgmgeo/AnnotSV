# AnnotSV: An integrated tool for Structural Variations annotation and ranking 

## QUICK INSTALLATION

1. The sources can be cloned to any directory:
```
cd /path/to/install/
git clone https://github.com/lgmgeo/AnnotSV.git
```
2. Then, the user can easily install the package using make:
```
cd /path/to/install/AnnotSV
make PREFIX=. install
make PREFIX=. install-human-annotation
make PREFIX=. install-mouse-annotation
```

3. Set the global environmental variable as the location of the git repo on your system. 

In csh:
```
setenv ANNOTSV /path/to/install/AnnotSV
```
In bash:
```
export ANNOTSV=/path/to/install/AnnotSV
```

## TEST

1. Change to the repo directory, and run the test
```
cd /path/to/install/AnnotSV/share/doc/AnnotSV/Example/
$ANNOTSV/bin/AnnotSV -SVinputFile test.bed -outputFile ./test.annotated.tsv -svtBEDcol 4
```
2. Examine the output

Happy exploring!


## COLLABORATIVE WORK

Anyone interested in implementing new annotations/features in AnnotSV?

Thanks to the [AnnotSV user community](https://lbgi.fr/AnnotSV/acknowledgments):

    - Bugs could be tackled efficiently

    - New ideas could be investigated faster


I look forward to the opportunity to work together,

feel free to fork the page if you want to help :-)

Véronique

