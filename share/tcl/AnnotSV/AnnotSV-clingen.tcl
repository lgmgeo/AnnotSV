############################################################################################################
# AnnotSV 3.4.6                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-present Veronique Geoffroy (veronique.geoffroy@inserm.fr)                             #
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


# ClinGenFileDownloaded:
########################
# ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv
# Some genes are on chrom X AND chrom Y:
# CSF2RA  1438    Xp22.33 and Yp11.2      chrX:1387693-1428828    30      Gene associated with autosomal recessive phenotype       0       No evidence available       2021-08-23       MONDO:0012580
# CSF2RA  1438    Xp22.33 and Yp11.2      chrY:1268814-1325218    30      Gene associated with autosomal recessive phenotype       0       No evidence available       2021-08-23       MONDO:0012580
# => be carreful to avoid redundancy


## - Check if the "ClinGen_gene_curation_list_GRCh37.tsv" file exists.
## - Check and create if necessary the 'date'_ClinGenAnnotations.tsv file.
proc checkClinGenFile {} {
    
    global g_AnnotSV
    
    
    set clingenDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/ClinGen"
    
    ## Check if the ClinGen file has been downloaded the formatted
    ##############################################################
    set ClinGenFileDownloaded [glob -nocomplain "$clingenDir/ClinGen_gene_curation_list_*.tsv"]
    set ClinGenFileFormattedGzip [glob -nocomplain "$clingenDir/*_ClinGenAnnotations.tsv.gz"]
    
    if {$ClinGenFileDownloaded eq "" && $ClinGenFileFormattedGzip eq ""} {
        # No ClinGene annotation
        return
    }
    
    if {[llength $ClinGenFileFormattedGzip]>1} {
        puts "Several ClinGen files exist:"
        puts "$ClinGenFileFormattedGzip"
        puts "Keep only one: [lindex $ClinGenFileFormattedGzip end]\n"
        foreach ClinGenF [lrange $ClinGenFileFormattedGzip 0 end-1] {
            file rename -force $ClinGenF $ClinGenF.notused
        }
        return
    }
    if {$ClinGenFileFormattedGzip eq ""} {
        ## - Create the 'date'_ClinGenAnnotations.tsv file.
        ##   Header: genes, HI and TS
        
        set ClinGenFileFormatted "$clingenDir/[clock format [clock seconds] -format "%Y%m%d"]_ClinGenAnnotations.tsv"
        puts "...creation of $ClinGenFileFormatted.gz ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])\n"
        ReplaceTextInFile "genes\tHI\tTS" $ClinGenFileFormatted
        
        # Parsing of $ClinGenFileDownloaded
        foreach L [LinesFromFile $ClinGenFileDownloaded] {
            set Ls [split $L "\t"]
            
            if {[regexp "^#Gene Symbol" $L]} {
                set i_gene [lsearch -exact $Ls "#Gene Symbol"];             if {$i_gene == -1} {puts "Bad header line syntax. Gene Symbol column not found - Exit with error"; exit 2}
                set i_HI   [lsearch -exact $Ls "Haploinsufficiency Score"]; if {$i_HI == -1} {puts "Bad header line syntax. Haploinsufficiency Score column not found - Exit with error"; exit 2}
                set i_TS   [lsearch -exact $Ls "Triplosensitivity Score"];  if {$i_TS == -1} {puts "Bad header line syntax. Triplosensitivity Score column not found - Exit with error"; exit 2}
                continue
            }
            if {[regexp "^#" $L]} {continue}
            
            set gene [lindex $Ls $i_gene]
            set HI [lindex $Ls $i_HI]
            set TS [lindex $Ls $i_TS]

			# To avoir redundancy for genes (CSF2RA and VAMP7) presents on chrom X AND chrom Y
			if {[info exists seen($gene)]} {continue}
			set seen($gene) 1

            if {$HI eq "Not yet evaluated"} {set HI ""}
            if {$TS eq "Not yet evaluated"} {set TS ""}
            
            lappend L_Texte "$gene\t$HI\t$TS"
        }
        
        # creation of $ClinGenFileFormatted.gz
        WriteTextInFile "[join $L_Texte "\n"]" $ClinGenFileFormatted
        if {[catch {exec gzip $ClinGenFileFormatted} Message]} {
            puts "-- checkClinGenFile --"
            puts "gzip $ClinGenFileFormatted"
            puts "$Message\n"
        }
        
        file delete -force $ClinGenFileDownloaded
    }
}

