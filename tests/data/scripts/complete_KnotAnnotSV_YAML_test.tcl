#!/usr/bin/env tclsh

############################################################################################################
# AnnotSV web site                                                                                         #
#                                                                                                          #
# Auteur: Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                                #
# Date: 2020/12/16                                                                                         #
############################################################################################################

proc WriteTextInFile {texte fichier} {
    set    fifi [open $fichier a]
    puts  $fifi $texte
    close $fifi
    return 1
}


# "./output/test.41_SV.annotated.tsv"
set AnnotSVoutput [lindex $argv 0]

# "/maison/geoffroy/knotAnnotSV/lgmgeoFork/knotAnnotSV/config_AnnotSV.yaml"
set yamlFile [lindex $argv 1]

# "./output/config_AnnotSV.tmp.yaml"
set yamlFileCompleted [lindex $argv 2]
file delete -force $yamlFileCompleted


set L_linesToWrite {}

# Look at the last POSITION used in the yamlFile
set lastPos 0
set f [open $yamlFile]
while {![eof $f]} {
    set L [gets $f]
    lappend L_linesToWrite $L
    if {[regexp "^ +POSITION: +(\[0-9\]+)" $L match pos]} {
	if {$pos > $lastPos} {set lastPos $pos}
    }
}

# Look at the user col name to add in "yamlFileCompleted"
set L_userColName {}
set f [open $AnnotSVoutput]
while {![eof $f]} {
    set L [gets $f]
    if {[regexp "^AnnotSV_ID" $L]} {
	if {![regsub "^.*SV_type\t" $L "" L]} {
	    regsub "^.*SV_length\t" $L "" L
	}
	regsub "Annotation_mode.*" $L "" L
	regsub "\t$" $L "" L
	set L_userColName [split $L "\t"]
	break
    }
}

# Complete the yamlFile (creation of yamlFileCompleted)
foreach col $L_userColName {
    incr lastPos
    lappend L_linesToWrite "${col}: # additional user column"
    lappend L_linesToWrite "    POSITION: $lastPos"
    lappend L_linesToWrite ""
}
WriteTextInFile [join $L_linesToWrite "\n"] $yamlFileCompleted


