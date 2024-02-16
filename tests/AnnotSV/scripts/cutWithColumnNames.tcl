#!/usr/bin/env tclsh

source /home/geoffroy/Tcl/general.tcl


## Command line example:
## /maison/geoffroy/Tcl/cutWithColumnNames.tcl test.annotated.tsv "SV type;Samples_ID;AnnotSV ranking"

#############
### INPUT ###
#############
# Fichier d'output d'AnnotSV :
set annotatedSVfileToParse [lindex $argv 0]
# Liste du nom des colonnes à afficher, séparées par des ";" (ET sans espace !)
set columnNames [lindex $argv 1]
set L_colNames [split $columnNames ";"]


# 1 - Cherche les numéros des colonnes à afficher 
# 2 - Affiche les infos dans les colonnes spécifiées en input
##############################################################
foreach L [LinesFromFile $annotatedSVfileToParse] {
	set Ls [split $L "\t"]
	if {[regexp "^## " $L]} {continue}

	# header
	if {[regexp "^AnnotSV_ID|^AnnotSV ID|^#|^gene|^variantID" $L]} {
		set lengthHeader [llength $Ls]
		set L_i {}
		set L_colOk {}
		foreach col $L_colNames {
			set i [lsearch $Ls "$col"] 
			if {$i ne -1} {	
				lappend L_i $i
				lappend L_colOk $col
			}
		}
		if {$L_colOk eq ""} {puts "No column to display";exit}
		puts "[join $L_colOk "\t"]"
		continue
	}

	# annotations
	if {[llength $Ls] ne $lengthHeader} {
		puts "different lengths of lines ([llength $L] ne $lengthHeader)"
	}
	foreach i $L_i {
		puts -nonewline "[lindex $Ls $i]\t"
	}
	puts ""

}


