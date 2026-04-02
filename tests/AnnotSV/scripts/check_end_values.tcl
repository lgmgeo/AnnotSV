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
    if {[regexp $svtype "TRA BND"]} {
        if {$start ne $end} {
            #puts "--- TRA/BND: $L"
            set test "TRAorBNDnotOK"
        }
    } else {
        if {$END ne $end} {
            puts "--- $END ne $end : $L"
            set test "SVnotOK"
        }
    }
}

puts $test
return $test

