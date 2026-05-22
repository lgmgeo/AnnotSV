############################################################################################################
# AnnotSV 3.5.8                                                                                            #
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

# SRC: Distribution directory (source)
# INSTALL: Installation directory (target)


SHELL = /usr/bin/env bash


DESTDIR               ?=
PREFIX                ?= /usr/local
INSTALLDIR1           := $(shell readlink -f "$(DESTDIR)$(PREFIX)")
INSTALLDIR2           := $(shell readlink -f "$(DESTDIR).")
BIN_INSTALL_DIR       := $(PREFIX)/bin
ETC_INSTALL_DIR       := $(PREFIX)/etc
SHARE_INSTALL_DIR     := $(PREFIX)/share
DOC_INSTALL_DIR       := $(SHARE_INSTALL_DIR)/doc
BASH_INSTALL_DIR      := $(SHARE_INSTALL_DIR)/bash
TESTS_INSTALL_DIR     := $(PREFIX)/tests
TCL_VERSION           := tcl$(shell echo 'puts $${tcl_version};exit 0' | tclsh)
TCL_SRC_DIR           := share/tcl
TCL_INSTALL_DIR       := $(SHARE_INSTALL_DIR)/$(TCL_VERSION)
PYTHON_INSTALL_DIR    := $(SHARE_INSTALL_DIR)/python3
JAR_INSTALL_DIR       := $(SHARE_INSTALL_DIR)/AnnotSV/jar
ANNOTSV_VERSION       := 3.5.8
HUMAN_VERSION         := 3.5
EXOMISER_VERSION      := 2406
MOUSE_VERSION         := 3.4.2
RM                    := /bin/rm
RMDIR                 := /bin/rmdir
MKDIR                 := install -d
MV                    := /bin/mv
CP                    := install -p -m 0644
CPDIR                 := /bin/cp -r
CHMOD                 := /bin/chmod -R 777
CONFIGFILE            := etc/AnnotSV/configfile
MAKEFILE              := Makefile
SRC_BASH_SCRIPTS      := $(shell find share/bash/AnnotSV/ -name '*.sh' 2> /dev/null)
SRC_DOCUMENTATIONS    := $(shell find License.txt changeLog.txt commandLineOptions.txt README.AnnotSV_*.pdf 2> /dev/null)
VC_INSTALL_FLAG       := $(DESTDIR)$(PYTHON_INSTALL_DIR)/variantconvert/pipinstall.flag
VC_VERSION            := 2.0.1
VC_CONFIG_INSTALL_DIR := $(DESTDIR)$(PYTHON_INSTALL_DIR)/variantconvert/src/variantconvert/configs
USE_ANNOTATIONS_DIR   := #flag whether separate annotation resources directory needed (e.g. for HPC environvment)
EX_SRC_PROPERTIES     := etc/AnnotSV/application.properties
EX_RP_FILE            := #optional filepath for previously downloaded EXomiser Rest-Prioritise

# make install
.PHONY: install

# make all
# make install                          (install-exomiser-resources is executed as part of install)
# make install install-human-annotation (install-exomiser is executed as part of install-human-annotation)
# make install install-mouse-annotation
ifeq ('$(INSTALLDIR1)' , '$(INSTALLDIR2)')
all: install-display install-documentationlight install-variantconvert install-exomiser-resources install-done
install: install-display install-documentationlight install-variantconvert install-exomiser-resources install-done
else
all: install-display install-configfile install-makefile install-executable install-tcl-toolbox install-bash-toolbox install-doc install-others-doc install-variantconvert install-exomiser-resources install-done
install: install-display install-configfile install-makefile install-executable install-tcl-toolbox install-bash-toolbox install-doc install-others-doc install-variantconvert install-exomiser-resources install-done
endif

# For all PREFIX
install-display:
	@echo ""
	@echo "Installation of AnnotSV-$(ANNOTSV_VERSION):"
	@echo "----------------------------"
	@echo DESTDIR=$(DESTDIR)
	@echo PREFIX=$(PREFIX)
	@echo TCL_VERSION=$(TCL_VERSION)

# For PREFIX==.
install-documentationlight: $(SRC_DOCUMENTATIONS)
	@echo ""
	@echo "Documentation light installation"
	@echo "--------------------------------"
	@if [ -n "$^" ]; then \
	    $(MV) $^ $(DESTDIR)$(DOC_INSTALL_DIR)/AnnotSV; \
	fi
	@if [ -d "$(TCL_SRC_DIR)" ]; then \
		$(MV) "$(TCL_SRC_DIR)" "$(TCL_INSTALL_DIR)"; \
	fi


# For PREFIX!=.
install-configfile: $(CONFIGFILE)
	@echo ""
	@echo "Configfile configuration"
	@echo "------------------------"
	$(MKDIR) $(DESTDIR)$(ETC_INSTALL_DIR)/AnnotSV
	install -p -m 0755 $(CONFIGFILE)  $(DESTDIR)$(ETC_INSTALL_DIR)/AnnotSV

# For PREFIX!=.
install-makefile: $(MAKEFILE)
	@echo ""
	@echo "Makefile installation"
	@echo "---------------------"
	install -p -m 0755 $(MAKEFILE)  $(DESTDIR)$(PREFIX)

# For PREFIX!=.
install-executable:
	@echo ""
	@echo "Executable installation"
	@echo "-----------------------"
	$(MKDIR) $(DESTDIR)$(BIN_INSTALL_DIR)
	install -p -m 0755 bin/AnnotSV $(DESTDIR)$(BIN_INSTALL_DIR)

# For PREFIX!=.
install-tcl-toolbox: 
	@echo ""
	@echo "Tcl scripts installation"
	@echo "------------------------"
	$(MKDIR) $(DESTDIR)$(TCL_INSTALL_DIR)/AnnotSV
	cd share/tcl ; tar cf - AnnotSV | tar xf - -C $(DESTDIR)$(TCL_INSTALL_DIR)/

# For all PREFIX
install-variantconvert:
	@echo ""
	@echo "variantconvert installation"
	@echo "---------------------------"
	
	@if [ -d $(DESTDIR)$(PYTHON_INSTALL_DIR)/variantconvert ]; then \
		echo "variantconvert directory found; purging locally before re-installing."; \
		rm -rf $(DESTDIR)$(PYTHON_INSTALL_DIR)/variantconvert/; \
	fi
	git clone --branch $(VC_VERSION) https://github.com/SamuelNicaise/variantconvert.git $(DESTDIR)$(PYTHON_INSTALL_DIR)/variantconvert/

	touch $(VC_INSTALL_FLAG)
	chmod 777 $(VC_INSTALL_FLAG)
	pip3 install -e $(DESTDIR)$(PYTHON_INSTALL_DIR)/variantconvert/. > ./tmp.variantconvert.txt 2>&1 \
	|| pip install -e $(DESTDIR)$(PYTHON_INSTALL_DIR)/variantconvert/. >> ./tmp.variantconvert.txt 2>&1 \
	|| python -m pip install -e $(DESTDIR)$(PYTHON_INSTALL_DIR)/variantconvert/. >> ./tmp.variantconvert.txt 2>&1 \
	|| rm -f $(VC_INSTALL_FLAG)
	@if [ -f $(VC_INSTALL_FLAG) ]; then \
		echo "variantconvert installed"; \
		$(CHMOD) ./tmp.variantconvert.txt; \
		rm -f ./tmp.variantconvert.txt; \
		$(CHMOD) $(VC_CONFIG_INSTALL_DIR); \
		$(MV) $(VC_CONFIG_INSTALL_DIR)/hs1 $(VC_CONFIG_INSTALL_DIR)/CHM13; \
		for f in $(VC_CONFIG_INSTALL_DIR)/CHM13/annotsv*; do \
			case "$$f" in \
				*.json) \
					sed -i 's/"##contig=<ID=chr/"##contig=<ID=/g' "$$f" ;; \
			esac; \
		done; \
		# Creation of the "*.local.json" files for Conda use. \
		for genome in GRCh37 GRCh38 CHM13; do \
			for source in bed vcf; do \
				for type in combined full fullsplit; do \
					echo "touch $(VC_CONFIG_INSTALL_DIR)/$$genome/annotsv3_from_$$source.$$type.local.json" ; \
					touch $(VC_CONFIG_INSTALL_DIR)/$$genome/annotsv3_from_$$source.$$type.local.json; \
					$(CHMOD) $(VC_CONFIG_INSTALL_DIR)/$$genome/annotsv3_from_$$source.$$type.local.json; \
				done; \
			done; \
		done; \
	else \
		echo "variantconvert not installed"; \
	fi


# For PREFIX!=.
install-bash-toolbox: $(SRC_BASH_SCRIPTS)
	@echo ""
	@echo "Bash scripts installation"
	@echo "-------------------------"
	$(MKDIR) $(DESTDIR)$(BASH_INSTALL_DIR)/AnnotSV
	$(CP) $^ $(DESTDIR)$(BASH_INSTALL_DIR)/AnnotSV

# For PREFIX!=.
install-doc: $(SRC_DOCUMENTATIONS)
	@echo ""
	@echo "Documentations installation"
	@echo "---------------------------"
	$(MKDIR) $(DESTDIR)$(DOC_INSTALL_DIR)/AnnotSV
	$(CP) $^ $(DESTDIR)$(DOC_INSTALL_DIR)/AnnotSV

# For PREFIX!=.
install-others-doc: share/doc/AnnotSV/Example
	$(CPDIR) $^ $(DESTDIR)$(DOC_INSTALL_DIR)/AnnotSV


# For all PREFIX
install-done: 
	@echo ""
	@echo "Done"
	@echo ""
	@echo "WARNING: Annotations need to be installed:"
	@echo "make DESTDIR=$(DESTDIR) PREFIX=$(PREFIX) install-human-annotation"
	@echo "make DESTDIR=$(DESTDIR) PREFIX=$(PREFIX) install-mouse-annotation"



# make install_organism_annotations
install-all-annotations: install-human-annotation install-mouse-annotation                                     

install-human-annotation: install-exomiser $(if $(USE_ANNOTATIONS_DIR),,Annotations_Human_$(HUMAN_VERSION).tar.gz)
ifndef USE_ANNOTATIONS_DIR
	@echo ""
	@echo "Installation of human annotation:"
	@echo ""
	tar -xf Annotations_Human_$(HUMAN_VERSION).tar.gz -C $(DESTDIR)$(SHARE_INSTALL_DIR)/AnnotSV/
	$(RM) -rf Annotations_Human_$(HUMAN_VERSION).tar.gz
	$(CHMOD) $(DESTDIR)$(SHARE_INSTALL_DIR)/AnnotSV/Annotations_*
	@echo ""
	@echo "--> Human annotation installed"
else
	@echo ""
	@echo "Flag for custom annotationDir; skipping local install of human annotations"
	@echo ""
endif


install-exomiser: $(if $(USE_ANNOTATIONS_DIR),,$(EXOMISER_VERSION)_phenotype.zip)
	@echo ""
	@echo "Installation of Exomiser data:"
	@echo ""
	
ifndef USE_ANNOTATIONS_DIR
	$(MKDIR) -p $(DESTDIR)$(SHARE_INSTALL_DIR)/AnnotSV/Annotations_Exomiser/$(EXOMISER_VERSION)
	unzip $(EXOMISER_VERSION)_phenotype.zip -d $(DESTDIR)$(SHARE_INSTALL_DIR)/AnnotSV/Annotations_Exomiser/$(EXOMISER_VERSION)/
	$(RM) -rf $(EXOMISER_VERSION)_phenotype.zip
else
	@echo ""
	@echo "Flag for custom annotationDir; skipped Exomiser phenotypes local installation"
	@echo ""
endif
	
# For all PREFIX
# application.properties
# exomiser-rest-prioritiser-14.1.0.jar
install-exomiser-resources: $(EX_SRC_PROPERTIES)
ifneq ('$(INSTALLDIR1)' , '$(INSTALLDIR2)')
	install -D -p -m 0755 $(EX_SRC_PROPERTIES) $(DESTDIR)$(ETC_INSTALL_DIR)/AnnotSV
endif

	$(MKDIR) -p $(DESTDIR)$(JAR_INSTALL_DIR)
ifndef EX_RP_FILE 
ifeq ($(wildcard $(DESTDIR)$(JAR_INSTALL_DIR)/exomiser-rest-prioritiser-14.1.0.jar),)
	curl -C - -LO https://github.com/exomiser/Exomiser/releases/download/14.1.0/exomiser-rest-prioritiser-14.1.0.jar
	install -p -m 0755 exomiser-rest-prioritiser-14.1.0.jar $(DESTDIR)$(JAR_INSTALL_DIR)/
	$(RM) exomiser-rest-prioritiser-14.1.0.jar
endif
else
	@echo "Custom Exomiser rest-priotiser path provided; creating symlink"
	ln -sf $(EX_RP_FILE) $(DESTDIR)$(JAR_INSTALL_DIR)/$(notdir $(EXRP_FILE))
endif

install-mouse-annotation: $(if $(USE_ANNOTATIONS_DIR),,Annotations_Mouse_$(MOUSE_VERSION).tar.gz) 
ifndef USE_ANNOTATIONS_DIR
	@echo ""
	@echo "Installation of mouse annotation:"
	@echo ""	
	$(MKDIR) $(DESTDIR)$(SHARE_INSTALL_DIR)/AnnotSV/
	tar -xf Annotations_Mouse_$(MOUSE_VERSION).tar.gz -C $(DESTDIR)$(SHARE_INSTALL_DIR)/AnnotSV/
	$(RM) -rf Annotations_Mouse_$(MOUSE_VERSION).tar.gz
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
	@echo "Uninstalling of AnnotSV"
	@echo "------------------------"
	$(RM) -f $(DESTDIR)$(BIN_INSTALL_DIR)/AnnotSV
	$(RM) -f $(DESTDIR)$(BIN_INSTALL_DIR)/INSTALL_annotations.sh
	$(RM) -f $(DESTDIR)$(BIN_INSTALL_DIR)/INSTALL_code.sh
	$(RM) -rf $(DESTDIR)$(TCL_INSTALL_DIR)/AnnotSV
	$(RM) -rf $(DESTDIR)$(PYTHON_INSTALL_DIR)/AnnotSV
	$(RM) -rf $(DESTDIR)$(PYTHON_INSTALL_DIR)/variantconvert
	$(RM) -rf $(DESTDIR)$(DOC_INSTALL_DIR)/AnnotSV
	$(RM) -rf $(DESTDIR)$(SHARE_INSTALL_DIR)/AnnotSV
	$(RM) -rf $(DESTDIR)$(BASH_INSTALL_DIR)/AnnotSV
	$(RM) -rf $(DESTDIR)$(ETC_INSTALL_DIR)/AnnotSV
	$(RM) -rf $(DESTDIR)$(TESTS_INSTALL_DIR)/AnnotSV
	$(RM) -rf $(DESTDIR)$(PREFIX)/Makefile
	$(RM) -rf $(DESTDIR)$(PREFIX)/README.md
	$(RM) -rf $(DESTDIR)$(PREFIX)/Scoring_Criteria_AnnotSV_*.xlsx
	$(RM) -rf $(DESTDIR)$(PREFIX)/.git
	$(RM) -rf $(DESTDIR)$(PREFIX)/.gitignore
	$(RM) -rf $(DESTDIR)$(PREFIX)/tmp.variantconvert.txt

uninstall2:
	$(RMDIR) --ignore-fail-on-non-empty $(DESTDIR)$(BIN_INSTALL_DIR) $(DESTDIR)$(BASH_INSTALL_DIR) $(DESTDIR)$(TCL_INSTALL_DIR) $(DESTDIR)$(PYTHON_INSTALL_DIR) $(DESTDIR)$(DOC_INSTALL_DIR) $(DESTDIR)$(SHARE_INSTALL_DIR) $(DESTDIR)$(ETC_INSTALL_DIR) $(DESTDIR)$(TESTS_INSTALL_DIR)

uninstall3:
	$(RMDIR) --ignore-fail-on-non-empty $(DESTDIR)$(PREFIX)

uninstall4:
	@echo ""
	@echo "Done"

