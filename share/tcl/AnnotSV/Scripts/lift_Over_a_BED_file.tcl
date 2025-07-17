#!/usr/bin/env tclsh


##########################################################################################################
# When performing a liftover on a BED file with more than 3 columns (with annotation columns), it may fail
# if some annotation values are empty.
# This script allows:
# - performing a liftover on the first three columns (chrom, start, end)
# - keeping associated annotation columns (even if empty)
# - adding the source of the chain file (e.g., "b38-GCA_009914755.4.over.chain" => source = "b38")
#   in front of annotations that represent genomic positions (e.g. 3:12345-45678 => b38_3:12345-45678)
##########################################################################################################

proc ContentFromFile {{Fichier ""}} {
    if {[string equal $Fichier ""]} {return ""}
    set f     [open $Fichier r]
    set Texte [read -nonewline $f]
    close $f
    return $Texte
}

proc LinesFromFile {{Fichier ""}} {
    return [split [ContentFromFile $Fichier] "\n"]
}

proc WriteTextInFile {texte fichier} {
    set    fifi [open $fichier a]
    puts  $fifi $texte
    close $fifi
    return 1
}


# Inputs
set BEDfileToLift [lindex $argv 0]
set outputFile [lindex $argv 1]
file delete -force $outputFile
set chainFile [lindex $argv 2]
set liftoverTool [lindex $argv 3]
if {$liftoverTool eq ""} {set liftoverTool "liftOver"}

if {[llength $argv] < 3} {
    puts "USAGE:\n$argv0 BEDfileToLift outputFile chainFile {liftoverTool}"
    exit
}


# Parsing the chain file name
# e.g:	b38-GCA_009914755.4.over.chain => "source = b38"
if {![regexp "^(.*?)(-|To|_)" [file tail $chainFile] match source]} {
    puts "...WARNING: No source extracted!"
} else {
    regsub -all "\"|'" $source "" source
    puts "...genome source: $source"
}


# Parsing of the BED file
set id 0
set L_toWrite {}
set pid [pid]
puts "...reading [file tail $BEDfileToLift]"
puts "...writing toLift.$pid.tmp.bed"
foreach L [LinesFromFile $BEDfileToLift] {
    if {[regexp "^#" $L]} {continue}
    if {$L == ""} {continue}
    set Ls [split $L "\t"]
    incr id
    
    # Creation of a BED file whith only the chrom, start and end values
    lappend L_toWrite "[join [lrange $Ls 0 2] "\t"]\t$id"
    if {![expr {$id%80000}]} {
        WriteTextInFile [join $L_toWrite "\n"] toLift.$pid.tmp.bed
        set L_toWrite {}
    }
    # Memorize the annotation columns
    set L_annotationsColumns {}
    foreach v [lrange $Ls 3 end] {
        if {[regexp "^\[0-9XYchrMT\]+:\[0-9\]+-\[0-9\]+" $v ]} {
            lappend L_annotationsColumns "${source}_$v"
        } else {
            lappend L_annotationsColumns "$v"
        }
    }
    set t($id) $L_annotationsColumns
}
WriteTextInFile [join $L_toWrite "\n"] toLift.$pid.tmp.bed


# Liftover
set command "$liftoverTool toLift.$pid.tmp.bed $chainFile lifted.$pid.tmp.bed unmapped.$pid.bed"
puts "...running:"
puts "$command"
if {[catch {eval exec $command} Message]} {
    puts "$Message"
    if {[regexp -nocase "error" $Message]} {
        file delete -force "./toLift.$pid.tmp.bed"
        file delete -force "./lifted.$pid.tmp.bed"
        file delete -force "./unmapped.$pid.bed"
        exit
    }
}


# Check and display different counts
set n_lifted [lindex [eval exec wc -l "./lifted.$pid.tmp.bed"] 0]
set n_unmapped [lindex [eval exec grep -c "^#" "./unmapped.$pid.bed"] 0]

puts "...$id genomic positions to lift"
puts "...$n_lifted genomic positions lifted"
puts "...$n_unmapped positions unmapped"

set a [expr {$id-$n_lifted}]
if {$a != $n_unmapped} {
    puts "\tplease check these counts"
}


# Writting the lifted output file
set L_toWrite {}
foreach L [LinesFromFile ./lifted.$pid.tmp.bed] {
    set Ls [split $L "\t"]
    set newLs [lrange $Ls 0 2]
    set id [lindex $Ls 3]
    lappend newLs {*}$t($id)
    lappend L_toWrite [join $newLs "\t"]
    if {![expr {$id%80000}]} {
        WriteTextInFile [join $L_toWrite "\n"] $outputFile.tmp
        set L_toWrite {}
    }
}
WriteTextInFile [join $L_toWrite "\n"] $outputFile.tmp


# Sort the lifted output file
eval exec sort -k1,1 -k2,2n $outputFile.tmp > $outputFile


# Clean
puts "...cleaning"
file delete -force "./toLift.$pid.tmp.bed"
file delete -force "./lifted.$pid.tmp.bed"
file delete -force "./unmapped.$pid.bed"
file delete -force $outputFile.tmp

