#!/usr/bin/env tclsh

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

set fileToCheck [lindex $argv 0]

set test "Ok"

foreach L [LinesFromFile $fileToCheck] {
        set INFO [lindex $L end]
        if {![regexp "END=(\[^;\]+)" $INFO match END]} {set END "---"}
        set svtype [lindex $L 1]
        set start [lindex $L end-2]
        set end [lindex $L end-1]
        if {![regexp "CIPOS(\[^=\])=(\[^;\]+)" [lindex $L end] match CIPOS]} {set CIPOS "---"}
        if {![regexp "CIEND(\[^=\])=(\[^;\]+)" [lindex $L end] match CIEND]} {set CIEND "---"}
        set cipos [lindex [split $CIPOS ","] 0]
        set ciend [lindex [split $CIEND ","] 1]
        if {$ciend eq ""} {continue}
        if {$svtype ne "TRA"} {
                if {[expr {$END+$ciend}] ne $end} {
                        puts "[lrange $L 0 4] -- cipos=$cipos ciend=$ciend END=$END"
			set test "notOK"
                }
        } else {
                # pour les TRA : end = $POS (du VCF) + $ciend + 1
                if {[expr {$end-$start}] ne [expr {$ciend-$cipos+1}]} {
                        puts "--- TRA: [lrange $L 0 4] -- [expr {$end-$start}] ne [expr {$ciend-$cipos+1}]"
			puts $L
                	set test "notOK"
		}
        }
}

puts $test
return $test



