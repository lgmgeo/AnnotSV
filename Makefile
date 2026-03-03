############################################################################################################
# AnnotSV 3.5.5                                                                                            #
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

SHELL = /usr/bin/env bash


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
TCLDIRDISTRIBUTED    := $(SHAREDIR)/tcl
TCLDIR               := $(SHAREDIR)/$(TCLVERSION)
PYTHONDIR            := $(SHAREDIR)/python3
ANNOTSV              := AnnotSV
JARDIR               := $(SHAREDIR)/$(ANNOTSV)/jar
VERSION              := 3.5.5
HUMANVERSION         := 3.5
MOUSEVERSION         := 3.4.2
RM                   := /bin/rm
RMDIR                := /bin/rmdir
MKDIR                := install -d
MV                   := /bin/mv
CP                   := install -p -m 0644
CPDIR                := /bin/cp -r
CHMOD                := /bin/chmod -R 777
CONFIGFILE           := $(ETCDIR)/$(ANNOTSV)/configfile
MAKEFILE             := Makefile
PROPERTIES           := $(ETCDIR)/$(ANNOTSV)/application.properties
BASH_SCRIPTS         := $(shell find $(BASHDIR)/$(ANNOTSV)/ -name '*.bash' 2> /dev/null)
DOCUMENTATIONS       := $(shell find License.txt changeLog.txt commandLineOptions.txt README.AnnotSV_*.pdf 2> /dev/null)
VC_FLAG              := $(DESTDIR)$(PYTHONDIR)/variantconvert/pipinstall.flag
VC_VERSION           := 2.0.1
VC_CONFIGDIR         := $(DESTDIR)$(PYTHONDIR)/variantconvert/src/variantconvert/configs
USEANNODIR           := #flag whether separate annotation resources directory needed (e.g. for HPC environvment)
EXRP_FILE            := #optional filepath for previously downloaded rest-prioritiser

# make install
.PHONY: install
ifeq ('$(INSTALLDIR1)' , '$(INSTALLDIR2)')
all: install-display install-documentationlight install-variantconvert install-done
install: install-display install-documentationlight install-variantconvert install-done
install-exomiser: install-exomiser-1 install-exomiser-3
else
all: install-display install-configfile install-makefile install-executable install-tcl-toolbox install-bash-toolbox install-doc install-others-doc install-variantconvert install-done
install: install-display install-configfile install-makefile install-executable install-tcl-toolbox install-bash-toolbox install-doc install-others-doc install-variantconvert install-done
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
	$(CP) $^ $(DESTDIR)$(DOCDIR)/$(ANNOTSV)
	$(CPDIR) $(TCLDIRDISTRIBUTED) $(TCLDIR)

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
	cd $(SHAREDIR)/tcl ; tar cf - $(ANNOTSV) | tar xf - -C $(DESTDIR)$(TCLDIR)/

install-variantconvert:
	@echo ""
	@echo "variantconvert installation"
	@echo "---------------------------"
	
	@if [ -d $(DESTDIR)$(PYTHONDIR)/variantconvert ]; then \
		echo "variantconvert directory found; purging locally before re-installing."; \
		rm -rf $(DESTDIR)$(PYTHONDIR)/variantconvert/; \
	fi
	git clone --branch $(VC_VERSION) https://github.com/SamuelNicaise/variantconvert.git $(DESTDIR)$(PYTHONDIR)/variantconvert/

	touch $(VC_FLAG)
	chmod 777 $(VC_FLAG)
	pip3 install -e $(DESTDIR)$(PYTHONDIR)/variantconvert/. > ./tmp.variantconvert.txt 2>&1 \
	|| pip install -e $(DESTDIR)$(PYTHONDIR)/variantconvert/. >> ./tmp.variantconvert.txt 2>&1 \
	|| python -m pip install -e $(DESTDIR)$(PYTHONDIR)/variantconvert/. >> ./tmp.variantconvert.txt 2>&1 \
	|| rm -f $(VC_FLAG)
	@if [ -f $(VC_FLAG) ]; then \
		echo "variantconvert installed"; \
		$(CHMOD) ./tmp.variantconvert.txt; \
		rm -f ./tmp.variantconvert.txt; \
		$(CHMOD) $(VC_CONFIGDIR); \
		$(MV) $(VC_CONFIGDIR)/hs1 $(VC_CONFIGDIR)/CHM13; \
		for f in $(VC_CONFIGDIR)/CHM13/annotsv*; do \
			case "$$f" in \
				*.json) \
					sed -i 's/"##contig=<ID=chr/"##contig=<ID=/g' "$$f" ;; \
			esac; \
		done; \
		# Creation of the "*.local.json" files for Conda use. \
		for genome in GRCh37 GRCh38 CHM13; do \
			for source in bed vcf; do \
				for type in combined full fullsplit; do \
					echo "touch $(VC_CONFIGDIR)/$$genome/annotsv3_from_$$source.$$type.local.json" ; \
					touch $(VC_CONFIGDIR)/$$genome/annotsv3_from_$$source.$$type.local.json; \
					$(CHMOD) $(VC_CONFIGDIR)/$$genome/annotsv3_from_$$source.$$type.local.json; \
				done; \
			done; \
		done; \
	else \
		echo "variantconvert not installed"; \
	fi

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

install-others-doc: $(DESTDIR)$(DOCDIR)/$(ANNOTSV)/Example
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

install-human-annotation: install-exomiser $(if $(USEANNODIR),,Annotations_Human_$(HUMANVERSION).tar.gz)
ifndef USEANNODIR
	@echo ""
	@echo "Installation of human annotation:"
	@echo ""
	tar -xf Annotations_Human_$(HUMANVERSION).tar.gz -C $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	$(RM) -rf Annotations_Human_$(HUMANVERSION).tar.gz
	$(CHMOD) $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_*
	@echo ""
	@echo "--> Human annotation installed"
else
	@echo ""
	@echo "Flag for custom annotationDir; skipping local install of human annotations"
	@echo ""
endif

install-exomiser-1: $(if $(USEANNODIR),,2406_phenotype.zip)
	@echo ""
	@echo "Installation of Exomiser data:"
	@echo ""
	
ifndef USEANNODIR
	$(MKDIR) -p $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_Exomiser/2406
	unzip 2406_phenotype.zip -d $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_Exomiser/2406/
	$(RM) -rf 2406_phenotype.zip
else
	@echo ""
	@echo "Flag for custom annotationDir; skipped Exomiser phenotypes local installation"
	@echo ""
endif
	
	$(MKDIR) -p $(DESTDIR)$(JARDIR)
ifndef EXRP_FILE
	curl -C - -LO https://github.com/exomiser/Exomiser/releases/download/14.1.0/exomiser-rest-prioritiser-14.1.0.jar
	install -p -m 0755 exomiser-rest-prioritiser-14.1.0.jar $(DESTDIR)$(JARDIR)/
	$(RM) exomiser-rest-prioritiser-14.1.0.jar
else
	@echo "Custom rest-priotiser path provided; creating symlink"
	ln -sf $(EXRP_FILE) $(DESTDIR)$(JARDIR)/$(notdir $(EXRP_FILE))
endif

install-exomiser-2:
	install -D -p -m 0755 $(PROPERTIES) $(DESTDIR)$(ETCDIR)/$(ANNOTSV)

install-exomiser-3:
	@echo ""
	@echo "--> Exomiser data installed"

install-mouse-annotation: $(if $(USEANNODIR),,Annotations_Mouse_$(MOUSEVERSION).tar.gz) 
ifndef USEANNODIR
	@echo ""
	@echo "Installation of mouse annotation:"
	@echo ""	
	$(MKDIR) $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	tar -xf Annotations_Mouse_$(MOUSEVERSION).tar.gz -C $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	$(RM) -rf Annotations_Mouse_$(MOUSEVERSION).tar.gz
	@echo ""
	@echo "--> Mouse annotation installed"
else
	@echo ""
	@echo "Flag for custom annotationDir; skipping local install of mouse annotations"
	@echo ""
endif


Annotations_%.tar.gz:
	@echo ""
	@echo "Download AnnotSV supporting data files:"
	@echo ""
	curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/$@

%_phenotype.zip:
	@echo ""
	@echo "Download Exomiser supporting data files:"
	@echo ""
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
	$(RM) -rf $(DESTDIR)$(PREFIX)/tmp.variantconvert.txt

uninstall2:
	$(RMDIR) --ignore-fail-on-non-empty $(DESTDIR)$(BINDIR) $(DESTDIR)$(BASHDIR) $(DESTDIR)$(TCLDIR) $(DESTDIR)$(PYTHONDIR) $(DESTDIR)$(DOCDIR) $(DESTDIR)$(SHAREDIR) $(DESTDIR)$(ETCDIR) $(DESTDIR)$(TESTSDIR)

uninstall3:
	$(RMDIR) --ignore-fail-on-non-empty $(DESTDIR)$(PREFIX)

uninstall4:
	@echo ""
	@echo "Done"

