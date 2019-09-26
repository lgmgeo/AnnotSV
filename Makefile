############################################################################################################
# AnnotSV 2.2.3                                                                                              #
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

SHELL = /bin/tcsh

DESTDIR              ?=
PREFIX               ?= /usr/local
BINDIR               ?= $(PREFIX)/bin
ETCDIR               ?= $(PREFIX)/etc
SHAREDIR             ?= $(PREFIX)/share
DOCDIR               ?= $(SHAREDIR)/doc
TCLVERSION            = tcl$(shell echo 'puts $${tcl_version};exit 0' | tclsh)
TCLDIRDISTRIBUTED     = share/tcl
TCLDIR               ?= $(SHAREDIR)/$(TCLVERSION)
ANNOTSV               = AnnotSV
VERSION               = 2.3
RM                    = /bin/rm
RMDIR                 = /bin/rmdir
MKDIR                 = install -d
MV                    = /bin/mv
CP                    = install -p -m 0644
CPDIR                 = /bin/cp -r
CONFIGFILE            = etc/$(ANNOTSV)/configfile
TCL_SCRIPTS           = $(shell find share/tcl/$(ANNOTSV)/ -name '*.tcl')
DOCUMENTATIONS        = License.txt changeLog.txt commandLineOptions.txt README.AnnotSV_$(VERSION).pdf
OTHERS_DOCUMENTATIONS = $(shell find share/doc/$(ANNOTSV)/ -type d -name 'Annotations_*')
ORGANISM              = Human


.PHONY: install
ifeq ('$(PREFIX)' , '.')
all: install-ligth
install: install-ligth
else
all: install-complete
install: install-complete
endif

install-ligth: $(DOCUMENTATIONS)
	$(MV) $^ $(DESTDIR)$(DOCDIR)/$(ANNOTSV)
	$(MV) $(TCLDIRDISTRIBUTED) $(TCLDIR)
	@echo "Done"

install-complete: install-configfile install-executable install-tcl-toolbox install-doc install-others-doc install-biological-data
	@echo ""
	@echo "Installation of $(ANNOTSV)-$(VERSION):"
	@echo "--------------------------------"
	@echo DESTDIR=$(DESTDIR)
	@echo PREFIX=$(PREFIX)
	@echo TCLVERSION=$(TCLVERSION)
	@echo ""
	@echo "Done"

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
	install -p -m 0755 bin/$(ANNOTSV)/AnnotSV.tcl $(DESTDIR)$(BINDIR)/$(ANNOTSV)

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

install-biological-data: share/doc/$(ANNOTSV)/Annotations_$(ORGANISM)
	$(MKDIR) $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/
	$(CPDIR) $^ $(DESTDIR)$(SHAREDIR)/$(ANNOTSV)/$(ORGANISM)


# make uninstall
.PHONY: uninstall
uninstall:
	@echo ""
	@echo "Uninstalling of $(ANNOTSV)"
	@echo "------------------------"
	$(RM) -r $(DESTDIR)$(BINDIR)/$(ANNOTSV)
	$(RM) -r $(DESTDIR)$(TCLDIR)/$(ANNOTSV)
	$(RM) -r $(DESTDIR)$(DOCDIR)/$(ANNOTSV)
	$(RM) -r $(DESTDIR)$(ETCDIR)/$(ANNOTSV)
	if ( $(PREFIX) != "/usr/local" ) then
		$(RMDIR) --ignore-fail-on-non-empty $(DESTDIR)$(BINDIR) $(DESTDIR)$(TCLDIR) $(DESTDIR)$(DOCDIR) $(DESTDIR)$(SHAREDIR) $(DESTDIR)$(ETCDIR) $(DESTDIR)$(PREFIX)
	endif
	@echo ""
	@echo "Done"



