############################################################################################################
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2019 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
#                                                                                                          #
# This program is free software; you can redistribute it and/or                                            #
# modify it under the terms of the GNU General Public License                                              #
# as published by the Free Software Foundation; either version 3                                           #
# of the License, or (at your option) any later version.                                                   #
#                                                                                                          #
# This program is distributed in the hope that it will be useful,                                          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                             #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################

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

