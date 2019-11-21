############################################################################################################
# AnnotSV 2.2.4                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2019 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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

SHELL = /bin/bash


DESTDIR              ?=
PREFIX               ?= /usr/local
INSTALLDIR1          := $(shell readlink -f "$(DESTDIR)$(PREFIX)")
INSTALLDIR2          := $(shell readlink -f "$(DESTDIR).")
BINDIR               := $(PREFIX)/bin
ETCDIR               := $(PREFIX)/etc
SHAREDIR             := $(PREFIX)/share
DOCDIR               := $(SHAREDIR)/doc
TCLVERSION           := tcl$(shell echo 'puts $${tcl_version};exit 0' | tclsh)
TCLDIRDISTRIBUTED    := share/tcl
TCLDIR               := $(SHAREDIR)/$(TCLVERSION)
ANNOTSV              := AnnotSV
VERSION              := 2.3
RM                   := /bin/rm
RMDIR                := /bin/rmdir
MKDIR                := install -d
MV                   := /bin/mv
CP                   := install -p -m 0644
CPDIR                := /bin/cp -r
CONFIGFILE           := etc/$(ANNOTSV)/configfile
PROPERTIES           := etc/$(ANNOTSV)/application.properties
TCL_SCRIPTS          := $(shell find share/tcl/$(ANNOTSV)/ -name '*.tcl')
DOCUMENTATIONS       := $(shell find License.txt changeLog.txt commandLineOptions.txt README.AnnotSV_*.pdf)


# make install
.PHONY: install
ifeq ('$(INSTALLDIR1)' , '$(INSTALLDIR2)')
all: install-display install-documentationlight install-done
install: install-display install-documentationlight install-done
else
all: install-display install-configfile install-executable install-tcl-toolbox install-doc install-others-doc install-done
install: install-display install-configfile install-executable install-tcl-toolbox install-doc install-others-doc install-done
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

install-executable:
	@echo ""
	@echo "Executable installation"
	@echo "-----------------------"
	$(MKDIR) $(DESTDIR)$(BINDIR)
	install -p -m 0755 bin/AnnotSV $(DESTDIR)$(BINDIR)

install-tcl-toolbox: $(TCL_SCRIPTS)
	@echo ""
	@echo "Tcl scripts installation"
	@echo "------------------------"
	$(MKDIR) $(DESTDIR)$(TCLDIR)/$(ANNOTSV)
	$(CP) $^ $(DESTDIR)$(TCLDIR)/$(ANNOTSV)

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
	tar xf Annotations_Human_$(VERSION).tar.gz -C $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	$(RM) -rf Annotations_Human_$(VERSION).tar.gz
	@echo ""
	@echo "--> Human annotation installed"

install-exomiser: 1902_phenotype.zip
	@echo ""
	@echo "Installation of Exomiser data:"
	@echo ""
	$(MKDIR) -p $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_Exomiser/1902
	install -p -m 0755 $(PROPERTIES) $(DESTDIR)$(ETCDIR)/$(ANNOTSV)
	$(CPDIR) share/AnnotSV/jar/ $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	unzip 1902_hg19.zip -d $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_Exomiser/1902/ -x 1902_hg19/1902_hg19_clinvar_whitelist.tsv.gz*
	unzip 1902_phenotype.zip -d $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/Annotations_Exomiser/1902/
	$(RM) -rf 1902_hg19.zip
	$(RM) -rf 1902_phenotype.zip
	@echo ""
	@echo "--> Exomiser data installed"

install-mouse-annotation: Annotations_Mouse_$(VERSION).tar.gz 
	@echo ""
	@echo "Installation of mouse annotation:"
	@echo ""
	$(MKDIR) $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	tar xf Annotations_Mouse_$(VERSION).tar.gz -C $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
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
	curl -C - -LO https://data.monarchinitiative.org/exomiser/data/1902_hg19.zip -LO https://data.monarchinitiative.org/exomiser/data/$@


# make uninstall
.PHONY: uninstall

ifeq ('$(PREFIX)' , '/usr/local')
uninstall: uninstall1 uninstall3
else
uninstall: uninstall1 uninstall2 uninstall3
endif

uninstall1:
	@echo ""
	@echo "Uninstalling of $(ANNOTSV)"
	@echo "------------------------"
	$(RM) -f $(DESTDIR)$(BINDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(TCLDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(DOCDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)
	$(RM) -rf $(DESTDIR)$(ETCDIR)/$(ANNOTSV)

uninstall2:
	$(RMDIR) --ignore-fail-on-non-empty $(DESTDIR)$(BINDIR) $(DESTDIR)$(TCLDIR) $(DESTDIR)$(DOCDIR) $(DESTDIR)$(SHAREDIR) $(DESTDIR)$(ETCDIR) $(DESTDIR)$(PREFIX)

uninstall3:
	@echo ""
	@echo "Done"

