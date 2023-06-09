<p align="center">
    <img src="docs/images/AnnotSV_logo.png" width="500">

<br />

# **AnnotSV: An integrated tool for Structural Variations annotation and ranking**

<br />

# Table of Contents

- ## [AnnotSV standalone README](/README.AnnotSV_3.3.4.pdf)
- ## [Webserver](docs/home.md)
- ## [Downloads](docs/downloads.md)
- ## [Ranking](docs/ranking.md)
- ## [Annotations](docs/annotations.md)

<br />

# Quick Installation

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

# TEST

1. Change to the repo directory, and run the test
```
cd /path/to/install/AnnotSV/share/doc/AnnotSV/Example/
$ANNOTSV/bin/AnnotSV -SVinputFile test.bed -outputFile ./test.annotated.tsv -svtBEDcol 4
```
2. Examine the output

Happy exploring!

<br />

# COLLABORATIVE WORK

Anyone interested in implementing new annotations/features in AnnotSV?

Thanks to the [AnnotSV user community](https://lbgi.fr/AnnotSV/acknowledgments):

    - Bugs could be tackled efficiently

    - New ideas could be investigated faster


I look forward to the opportunity to work together,

feel free to fork the page if you want to help :-)

Véronique

