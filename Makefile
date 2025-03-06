############################################################################################################
# AnnotSV 3.4.6                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-present Veronique Geoffroy (veronique.geoffroy@inserm.fr)                             #
#                                                                                                          #
# This is part of AnnotSV source code.                                                                     #
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

SHELL = /usr/bin/bash


DESTDIR              ?=
PREFIX               ?= /usr/local
INSTALLDIR1          := $(shell readlink -f "$(DESTDIR)$(PREFIX)")
INSTALLDIR2          := $(shell readlink -f "$(DESTDIR).")
BINDIR               := $(PREFIX)/bin
ETCDIR               := $(PREFIX)/etc
SHAREDIR             := $(PREFIX)/share
DOCDIR               := $(SHAREDIR)/doc
BASHDIR              := $(SHAREDIR)/bash
TESTSDIR             := $(PREFIX)/tests
TCLVERSION           := tcl$(shell echo 'puts $${tcl_version};exit 0' | tclsh)
TCLDIRDISTRIBUTED    := share/tcl
TCLDIR               := $(SHAREDIR)/$(TCLVERSION)
PYTHONDIR            := $(SHAREDIR)/python3
ANNOTSV              := AnnotSV
VERSION              := 3.4.6
RM                   := /bin/rm
RMDIR                := /bin/rmdir
MKDIR                := install -d
MV                   := /bin/mv
CP                   := install -p -m 0644
CPDIR                := /bin/cp -r
CHMOD                := /bin/chmod -R 777
CONFIGFILE           := etc/$(ANNOTSV)/configfile
MAKEFILE             := Makefile
PROPERTIES           := etc/$(ANNOTSV)/application.properties
BASH_SCRIPTS         := $(shell find share/bash/$(ANNOTSV)/ -name '*.bash' 2> /dev/null)
DOCUMENTATIONS       := $(shell find License.txt changeLog.txt commandLineOptions.txt README.AnnotSV_*.pdf 2> /dev/null)

# make install
.PHONY: install
ifeq ('$(INSTALLDIR1)' , '$(INSTALLDIR2)')
all: install-display install-documentationlight install-variantconvert install-done
install: install-display install-documentationlight install-variantconvert install-done
install-exomiser: install-exomiser-1 install-exomiser-3
else
all: install-display install-configfile install-makefile install-executable install-tcl-toolbox install-python-toolbox install-bash-toolbox install-doc install-others-doc install-variantconvert install-done
install: install-display install-configfile install-makefile install-executable install-tcl-toolbox install-python-toolbox install-bash-toolbox install-doc install-others-doc install-variantconvert install-done
install-exomiser: install-exomiser-1 install-exomiser-2 install-exomiser-3
endif

install-display:
	@echo ""
	@echo "Installation of $(ANNOTSV)-$(VERSION):"
	@echo "----------------------------"
	@echo DESTDIR=$(DESTDIR)
	@echo PREFIX=$(PREFIX)
	@echo TCLVERSION=$(TCLVERSION)

install-documentationlight: $(DOCUMENTATIONS)
	@echo ""
	$(MV) $^ $(DESTDIR)$(DOCDIR)/$(ANNOTSV)
	$(MV) $(TCLDIRDISTRIBUTED) $(TCLDIR)

install-configfile: $(CONFIGFILE)
	@echo ""
	@echo "Configfile configuration"
	@echo "------------------------"
	$(MKDIR) $(DESTDIR)$(ETCDIR)/$(ANNOTSV)
	install -p -m 0755 $(CONFIGFILE)  $(DESTDIR)$(ETCDIR)/$(ANNOTSV)

install-makefile: $(MAKEFILE)
	@echo ""
	@echo "Makefile installation"
	@echo "---------------------"
	install -p -m 0755 $(MAKEFILE)  $(DESTDIR)$(PREFIX)

install-executable:
	@echo ""
	@echo "Executable installation"
	@echo "-----------------------"
	$(MKDIR) $(DESTDIR)$(BINDIR)
	install -p -m 0755 bin/AnnotSV $(DESTDIR)$(BINDIR)

install-tcl-toolbox: 
	@echo ""
	@echo "Tcl scripts installation"
	@echo "------------------------"
	$(MKDIR) $(DESTDIR)$(TCLDIR)/$(ANNOTSV)
	cd share/tcl ; tar cf - $(ANNOTSV) | tar xf - -C $(DESTDIR)$(TCLDIR)/

install-python-toolbox:
	@echo ""
	@echo "Python scripts installation"
	@echo "---------------------------"
	$(MKDIR) $(DESTDIR)$(PYTHONDIR)/variantconvert
	cd share/python3 ; tar cf - variantconvert | tar xf - -C $(DESTDIR)$(PYTHONDIR)/

install-variantconvert:
	@echo ""
	@echo "variantconvert installation"
	@echo "---------------------------"
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/pipinstall.flag
	chmod 777 $(DESTDIR)$(PYTHONDIR)/variantconvert/pipinstall.flag
	pip3 install -e $(DESTDIR)$(PYTHONDIR)/variantconvert/. &> ./tmp.variantconvert.txt || pip install -e $(DESTDIR)$(PYTHONDIR)/variantconvert/. &> ./tmp.variantconvert.txt || rm -f $(DESTDIR)$(PYTHONDIR)/variantconvert/pipinstall.flag
	rm -f ./tmp.variantconvert.txt
	ls -l $(DESTDIR)$(PYTHONDIR)/variantconvert/pipinstall.flag || echo "variantconvert not installed"
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_bed.combined.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_bed.full.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_bed.fullsplit.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_vcf.combined.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_vcf.full.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_vcf.fullsplit.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_bed.combined.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_bed.full.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_bed.fullsplit.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_vcf.combined.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_vcf.full.local.json
	touch $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_vcf.fullsplit.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_bed.combined.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_bed.full.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_bed.fullsplit.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_vcf.combined.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_vcf.full.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh37/annotsv3_from_vcf.fullsplit.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_bed.combined.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_bed.full.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_bed.fullsplit.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_vcf.combined.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_vcf.full.local.json
	$(CHMOD) $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs/GRCh38/annotsv3_from_vcf.fullsplit.local.json


install-bash-toolbox: $(BASH_SCRIPTS)
	@echo ""
	@echo "Bash scripts installation"
	@echo "-------------------------"
	$(MKDIR) $(DESTDIR)$(BASHDIR)/$(ANNOTSV)
	$(CP) $^ $(DESTDIR)$(BASHDIR)/$(ANNOTSV)

install-doc: $(DOCUMENTATIONS)
	@echo ""
	@echo "Documentations installation"
	@echo "---------------------------"
	$(MKDIR) $(DESTDIR)$(DOCDIR)/$(ANNOTSV)
	$(CP) $^ $(DESTDIR)$(DOCDIR)/$(ANNOTSV)

install-others-doc: share/doc/$(ANNOTSV)/Example
	$(CPDIR) $^ $(DESTDIR)$(DOCDIR)/$(ANNOTSV)

install-done: 
	@echo ""
	@echo "Done"
	@echo ""
	@echo "WARNING: Annotations need to be installed:"
	@echo "make DESTDIR=$(DESTDIR) PREFIX=$(PREFIX) install-human-annotation"
	@echo "make DESTDIR=$(DESTDIR) PREFIX=$(PREFIX) install-mouse-annotation"



# make install_organism_annotations
install-all-annotations: install-human-annotation install-mouse-annotation                                     

install-human-annotation: Annotations_Human_$(VERSION).tar.gz install-exomiser
	@echo ""
	@echo "Installation of human annotation:"
	@echo ""
	tar -xf Annotations_Human_$(VERSION).tar.gz -C $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	$(RM) -rf Annotations_Human_$(VERSION).tar.gz
	$(CHMOD) $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_*
	@echo ""
	@echo "--> Human annotation installed"

install-exomiser-1: 2309_phenotype.zip
	@echo ""
	@echo "Installation of Exomiser data:"
	@echo ""
	$(MKDIR) -p $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_Exomiser/2309
	tar -xf 2309_hg19.tar.gz -C $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_Exomiser/2309/
	unzip 2309_phenotype.zip -d $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_Exomiser/2309/
	$(RM) -rf 2309_phenotype.zip
	$(RM) -rf 2309_hg19.tar.gz

install-exomiser-2:
	install -p -m 0755 $(PROPERTIES) $(DESTDIR)$(ETCDIR)/$(ANNOTSV)
	$(CPDIR) share/AnnotSV/jar/ $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/

install-exomiser-3:
	@echo ""
	@echo "--> Exomiser data installed"

install-mouse-annotation: Annotations_Mouse_$(VERSION).tar.gz 
	@echo ""
	@echo "Installation of mouse annotation:"
	@echo ""
	$(MKDIR) $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	tar -xf Annotations_Mouse_$(VERSION).tar.gz -C $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	$(RM) -rf Annotations_Mouse_$(VERSION).tar.gz
	@echo ""
	@echo "--> Mouse annotation installed"

Annotations_%.tar.gz:
	@echo ""
	@echo "Download AnnotSV supporting data files:"
	@echo ""
	curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/$@

%_phenotype.zip:
	@echo ""
	@echo "Download Exomiser supporting data files:"
	@echo ""
	curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/2309_hg19.tar.gz
	curl -C - -LO https://data.monarchinitiative.org/exomiser/data/$@


# make uninstall
.PHONY: uninstall

ifeq ('$(PREFIX)' , '/usr/local')
uninstall: uninstall1 uninstall4
else ifeq ('$(INSTALLDIR1)' , '$(INSTALLDIR2)')
uninstall: uninstall1 uninstall2 uninstall4
else
uninstall: uninstall1 uninstall2 uninstall3 uninstall4
endif

uninstall1:
	@echo ""
	@echo "Uninstalling of $(ANNOTSV)"
	@echo "------------------------"
	$(RM) -f $(DESTDIR)$(BINDIR)/$(ANNOTSV)
	$(RM) -f $(DESTDIR)$(BINDIR)/INSTALL_annotations.sh
	$(RM) -f $(DESTDIR)$(BINDIR)/INSTALL_code.sh
	$(RM) -rf $(DESTDIR)$(TCLDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(PYTHONDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(PYTHONDIR)/variantconvert
	$(RM) -rf $(DESTDIR)$(DOCDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(BASHDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(ETCDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(TESTSDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(PREFIX)/Makefile
	$(RM) -rf $(DESTDIR)$(PREFIX)/README.md
	$(RM) -rf $(DESTDIR)$(PREFIX)/Scoring_Criteria_AnnotSV_*.xlsx
	$(RM) -rf $(DESTDIR)$(PREFIX)/.git
	$(RM) -rf $(DESTDIR)$(PREFIX)/.gitignore

uninstall2:
	$(RMDIR) --ignore-fail-on-non-empty $(DESTDIR)$(BINDIR) $(DESTDIR)$(BASHDIR) $(DESTDIR)$(TCLDIR) $(DESTDIR)$(PYTHONDIR) $(DESTDIR)$(DOCDIR) $(DESTDIR)$(SHAREDIR) $(DESTDIR)$(ETCDIR) $(DESTDIR)$(TESTSDIR)

uninstall3:
	$(RMDIR) --ignore-fail-on-non-empty $(DESTDIR)$(PREFIX)

uninstall4:
	@echo ""
	@echo "Done"

