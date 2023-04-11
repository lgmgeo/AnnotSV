############################################################################################################
# AnnotSV 3.3.4                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2023 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                            #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################

proc isCandidate {gene} {
    global g_AnnotSV
    global g_candidate
    
    if {![info exists g_candidate(DONE)]} {
	set g_candidate(DONE) 1 
	
        # List of the candidate genes (given by the user)
	set L_Candidates ""
        if {$g_AnnotSV(candidateGenesFile) ne ""} {
	    set f [open $g_AnnotSV(candidateGenesFile)]
	    while {![eof $f]} {
	        set L [gets $f]
	        lappend L_Candidates {*}[split $L " |\n|\t"]
    	    }
	    close $f
        }
	foreach g $L_Candidates {
	    set g_candidate($g) 1
	}
    }
    
    if {[info exists g_candidate($gene)]} {
	return 1
    } else {
	return 0
    } 
}




