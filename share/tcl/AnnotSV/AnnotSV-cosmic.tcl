############################################################################################################
# AnnotSV 3.4                                                                                              #
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


## - Check if the "CosmicCompleteCNA.tsv.gz" file has been downloaded:
##   If yes, create the "CosmicCompleteCNA_$g_AnnotSV(genomeBuild).sorted.bed" file (in the "Users" directory) and delete the previous one.
proc checkCOSMICfile {} {
    
    global g_AnnotSV
    
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set COSMICdir "$extannDir/FtIncludedInSV/COSMIC"
    set COSMICfileDownloaded "$COSMICdir/$g_AnnotSV(genomeBuild)/CosmicCompleteCNA.tsv.gz"
    set COSMICfileFormatted "$extannDir/Users/$g_AnnotSV(genomeBuild)/FtIncludedInSV/CosmicCompleteCNA_$g_AnnotSV(genomeBuild).bed"
    
    if {[file exists $COSMICfileFormatted]} {
        set g_AnnotSV(COSMICann) 1
    } elseif {[file exists $COSMICfileDownloaded]} {
        # Creation of "CosmicCompleteCNA_$g_AnnotSV(genomeBuild).sorted.bed"
        puts "\t...COSMIC configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
        set i 0
        set f [open "| gzip -cd $COSMICfileDownloaded"]
        puts "\t\t...reading $COSMICfileDownloaded"
        while {![eof $f]} {
            set L [gets $f]
            set Ls [split $L "\t"]
            incr i
            if {[expr {$i%10000000}] eq 0} {puts "\t$i lines read"}
            if {[regexp "^COSMIC_CNV_ID|CNV_ID" $L]} {
                set i_id [lsearch -regexp $Ls "^COSMIC_CNV_ID|CNV_ID"]
                set i_type [lsearch -exact $Ls "MUT_TYPE"]
                set i_coord [lsearch -regexp $Ls "Chromosome:G_Start"] ;# Old raw data version (before 2023)
                set i_chrom [lsearch -regexp $Ls "CHROMOSOME"]    ;# New raw data version
                set i_start [lsearch -regexp $Ls "GENOME_START"]  ;# New raw data version
                set i_end [lsearch -regexp $Ls "GENOME_STOP"]     ;# New raw data version
                set L_texteToWrite(Header) "#Chrom\tstart\tend\tCOSMIC_ID\tCOSMIC_MUT_TYP"
                continue
            }
            if {$L eq ""} {continue}
            
            set id [lindex $Ls $i_id]
            set cnvtype [lindex $Ls $i_type]
            
            if {$i_coord ne -1} {
                set coord [lindex $Ls $i_coord]
                regexp "^(\[0-9\]+):(\[0-9\]+)\\.\\.(\[0-9\]+)$" $coord match chrom start end
            } else {
                set chrom [lindex $Ls $i_chrom]
                set start [lindex $Ls $i_start]
                set end   [lindex $Ls $i_end]
            }
            
            lappend L_texteToWrite($chrom) "$chrom\t$start\t$end\t$id\t$cnvtype"
        }
        
        puts "\t\t...creation of the COSMIC bed file in $extannDir/Users/$g_AnnotSV(genomeBuild)/FtIncludedInSV/"
        WriteTextInFile $L_texteToWrite(Header) $COSMICfileFormatted
        foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
            if {![info exists L_texteToWrite($chrom)]} {continue}
            
            ## Remove identical lines (due to liftover_to_hg19 for example)
            set L_texteToWrite($chrom) [lsort -unique $L_texteToWrite($chrom)]
            
            ## Put in karyotypic order
            set L_texteToWrite($chrom) [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_texteToWrite($chrom)]]
            
            ## Write the BED file
            WriteTextInFile [join $L_texteToWrite($chrom) "\n"] $COSMICfileFormatted
            unset L_texteToWrite($chrom)
        }
        
        file delete -force $COSMICfileDownloaded
        
        puts "\t\t   done ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
        set g_AnnotSV(COSMICann) 1
    }
    
    return
}

