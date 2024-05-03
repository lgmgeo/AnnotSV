############################################################################################################
# AnnotSV 3.4.1                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2024 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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


# - Check if the gnomAD (LOEUF) file has been downloaded:
#   "gnomad.v2.1.1.lof_metrics.by_gene.txt"
#
# - Check and create if necessary the following file:
#   'date'_gnomAD.LOEUF.pLI.annotations.tsv
proc checkLOEUFfile {} {
    
    global g_AnnotSV
    
    ## Check if the LOEUF file has been downloaded then formatted
    #############################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based"
    
    set LOEUFfileDownloaded [glob -nocomplain "$extannDir/gnomAD/gnomad.*.lof_metrics.by_gene.txt"]
    set LOEUFformattedFile  [glob -nocomplain "$extannDir/gnomAD/*_gnomAD.LOEUF.pLI.annotations.tsv.gz"]
    
    if {$LOEUFfileDownloaded eq "" && $LOEUFformattedFile eq ""} {
        # No "LOEUF" annotation
        return
    }
    
    if {[llength $LOEUFformattedFile]>1} {
        puts "Several LOEUF files exist:"
        puts "$LOEUFformattedFile"
        puts "Keep only one: [lindex $LOEUFformattedFile end]\n"
        foreach lf [lrange $LOEUFformattedFile 0 end-1] {
            file rename -force $lf $lf.notused
        }
        return
    }
    
    if {$LOEUFformattedFile eq ""} {
        ## Create : 'date'_gnomAD.LOEUF.pLI.annotations.tsv.gz   ; # Header: genes  LOEUF_bin
        set LOEUFformattedFile "$extannDir/gnomAD/[clock format [clock seconds] -format "%Y%m%d"]_gnomAD.LOEUF.pLI.annotations.tsv"
        
        puts "\t...creation of $LOEUFformattedFile.gz ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
        puts "\t   (done only once during the first LOEUF annotation)"
        
        ReplaceTextInFile "genes\tLOEUF_bin\tpLI_gnomAD\tpLI_ExAC" $LOEUFformattedFile
        foreach L [LinesFromFile $LOEUFfileDownloaded] {
            set Ls [split $L "\t"]
            if {[regexp  "^gene" $L]} {
                set i_gene    0
                set i_loeuf   [lsearch -regexp $Ls "^oe_lof_upper_bin"]; if {$i_loeuf == -1} {puts "Bad syntax into $LOEUFfileDownloaded.\noe_lof_upper_bin field not found - Exit with error"; exit 2}
                set i_pLI     [lsearch -regexp $Ls "^pLI"];              if {$i_pLI == -1}   {puts "Bad syntax into $LOEUFfileDownloaded.\npLI field not found - Exit with error"; exit 2}
                set i_pLIexac [lsearch -regexp $Ls "^exac_pLI"];         if {$i_pLIexac == -1} {puts "Bad syntax into $LOEUFfileDownloaded.\nexac_pLI field not found - Exit with error"; exit 2}
                continue
            }
            set gene    [lindex $Ls $i_gene]
            set loeuf   [lindex $Ls $i_loeuf]
            set pLI     [lindex $Ls $i_pLI]
            set pLIexac [lindex $Ls $i_pLIexac]
            if {$loeuf != "NA"}   {lappend L_loeuf($gene) "$loeuf"}     ;# Low LOEUF scores (e.g. 0 from values [0-9]) suggest a higher intolerance to inactivation
            if {$pLI != "NA"}     {lappend L_pLI($gene) "$pLI"}         ;# High pLI scores (e.g. 1 from values [0-1]) suggest a higher intolerance to inactivation
            if {$pLIexac != "NA"} {lappend L_pLIexac($gene) "$pLIexac"} ;# High pLIexac scores (e.g. 0) suggest a higher intolerance to inactivation
        }
        
        set TexteToWrite ""
        set allGenes [array names L_loeuf]
        lappend allGenes {*}[array names L_pLI]
        lappend allGenes {*}[array names L_pLIexac]
        set allGenes [lsort -unique $allGenes]
        foreach g $allGenes {
            if {![info exists L_loeuf($g)]} {set loeuf ""} else {set loeuf [lindex [lsort $L_loeuf($g)] 0]}
            if {![info exists L_pLI($g)]} {set pLI ""} else {set pLI [lindex [lsort $L_pLI($g)] 0]}
            if {![info exists L_pLIexac($g)]} {set pLIexac ""} else {set pLIexac [lindex [lsort $L_pLIexac($g)] 0]}
            
            lappend TexteToWrite "$g\t$loeuf\t$pLI\t$pLIexac"
        }
        WriteTextInFile [join $TexteToWrite "\n"] $LOEUFformattedFile
        if {[catch {exec gzip $LOEUFformattedFile} Message]} {
            puts "-- checkLOEUFfile --"
            puts "gzip $LOEUFformattedFile"
            puts "$Message\n"
        }
        
        file delete -force $LOEUFfileDownloaded
    }
}

