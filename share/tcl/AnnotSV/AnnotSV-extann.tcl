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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                            #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################


proc ExternalAnnotations args {

    # To add external annotation to a gene in particular
    # All .../Annotations_$g_AnnotSV(organism)/*/*.tsv and .../Annotations_$g_AnnotSV(organism)/*/*.tsv.gz files are used (defined in AnnotSV-main.tcl)
    #
    # Format is tab separated values, 1st line is a header, 1st "genes" column is the gene name, rest is free
    # Typical use would be a gene file containing specific annotations such as tranmission mode, disease, expression...
    # WARNING: shouldn't contain "{" and "}". Replaced here by "(" and ")"
    
    # args :
    #   - $F
    #   - "L_Files"    -return-> Annotation files list
    #   - $F,"Loaded"  -return-> 1
    #   - $F,"L_ID"    -return-> List of all genes from $F
    #   - $F,Header    -return-> Header from $F
    #   - $F,$ID       -return-> Annotation for $ID (= gene) from $F ([lrange $L 1 end])

    global g_AnnotSV
    global g_ExtAnnotation

    set What [join $args ","]
    if {[info exists g_ExtAnnotation($What)]} {return [set g_ExtAnnotation($What)]}

    if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) ne ""} {
	set L_Files_Anno [set g_AnnotSV(extann)]
    } else {return ""}
    
    if {[file exists [lindex $args 0]]} {
	set What [join $args ","]
	set L_Files_Anno [lindex $args 0]
    } else {
	set What [join [concat [list $L_Files_Anno] $args] ","]
    }
    
    #puts "--args $args"
    #puts "--$L_Files_Anno"
    #puts "--$What"

    if {[info exists g_ExtAnnotation($What)]} {return [set g_ExtAnnotation($What)]}
    foreach File_Anno [lsort -unique $L_Files_Anno] {
	if {![info exists g_ExtAnnotation($File_Anno,Loaded)]} {
	    
	    if {![info exists g_ExtAnnotation(L_Files)]} {
		set g_ExtAnnotation(L_Files) {}
	    }
	    lappend g_ExtAnnotation(L_Files) $File_Anno

	    set First 1
	    
	    set NbHeaderColumns 0
	    
	    set g_ExtAnnotation($File_Anno,Loaded) 1
	    
	    lappend g_ExtAnnotation(display) "\t\t...[file tail $File_Anno]"
	    
	    if {![info exists g_ExtAnnotation($File_Anno,L_ID)]} {
		set g_ExtAnnotation($File_Anno,L_ID) {}
	    }
	    if {[regexp ".gz$" $File_Anno]} {
		set lLines [LinesFromGZFile $File_Anno]
	    } else {
		set lLines [LinesFromFile $File_Anno]
	    }
	    foreach L $lLines {
		## Bug during WritingByGene if "{ or }" are presents.
		## To be done here, inside each element of the list (not on "$lLines").
		regsub -all "{" $L "(" L
		regsub -all "}" $L ")" L
		
		if {! [regexp -nocase {[a-z0-9]+} $L] || [regexp -nocase {^[\#]} $L]} {continue}
		
		if {$First} {
		    set First 0
		    if {![regexp -nocase "gene" [lindex [split $L "\t"] 0]]} {
			puts "\t\tReading $File_Anno:\n\tWARNING: Header is present ($L) but does not contain a \"gene\" named column."
			set g_ExtAnnotation($File_Anno,Header) ""
		    } else {
			set g_ExtAnnotation($File_Anno,Header) $L
			set NbHeaderColumns [llength [split $L "\t"]]
		    }
		    continue
		} 
		
		set Ls [split $L "\t"]
		set ID [string trim [lindex $Ls 0]]
		
		if {[llength $Ls] ne $NbHeaderColumns} {puts "WARNING: $ID, format (nb columns) is not good [llength $L] vs $NbHeaderColumns columns"}
		
		if {![info exists g_ExtAnnotation($File_Anno,$ID)]} {
		    lappend g_ExtAnnotation($File_Anno,L_ID) $ID
		} else {
		    puts "WARNING: Reading Annotation file [file tail $File_Anno] and $ID is seen multiple times."
		}
		
		set g_ExtAnnotation($File_Anno,$ID) [join [lrange $Ls 1 end] "\t"]
		
	    }
	    
	    lappend g_ExtAnnotation(display) "\t\t([llength [set g_ExtAnnotation($File_Anno,L_ID)]] gene identifiers and [expr {[llength [split [set g_ExtAnnotation($File_Anno,Header)] "\t"]]-1}] annotations columns: [join [lrange [split [set g_ExtAnnotation($File_Anno,Header)] "\t"] 1 end] ", "])"
	    
	    if {[info exists g_ExtAnnotation($What)]} {return [set g_ExtAnnotation($What)]}
	} 
    }
    return ""
}
