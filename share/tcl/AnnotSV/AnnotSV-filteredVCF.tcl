############################################################################################################
# AnnotSV 3.0.6                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2021 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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



##########################################################################
# Example with "-candidateSnvIndelSamples" = "sample1 sample2"
# -> Prepare 2 new annotation columns: "compound-htz($sample1) compound-htz($sample2)"
#    = List of htz variants (${chrom}_$pos) presents in the gene overlapped by the SV
##########################################################################


proc filteredVCFannotation {GENEchrom GENEstart GENEend Line headerOutput} {

    global g_AnnotSV
    global L_htzPos
    global L_samplesFromTheHeader

    ## VCFs parsing is done only 1 time (After, g_AnnotSV(filteredVCFparsing) is set to "done")
    if {![info exists g_AnnotSV(filteredVCFparsing)]} {
	set g_AnnotSV(filteredVCFparsing) "done"
	# parsing of $g_AnnotSV(candidateSnvIndelFiles): creation of L_htzPos($chrom,$sample)
	puts "\n\n...parsing of candidateSnvIndelFiles for \"$g_AnnotSV(candidateSnvIndelSamples)\" ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 
	
	# No htz compound for the full lines
	set L_htzPos(FULL) [lrepeat [llength $g_AnnotSV(candidateSnvIndelSamples)] ""]
	
	# "eval glob" accept regular expression ("*.vcf) as well as a list of files ("sample1.vcf sample2.vcf.gz"):
	foreach vcfF [eval glob -nocomplain $g_AnnotSV(candidateSnvIndelFiles)] {
	    puts "\t...parse all positions from $vcfF ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 
	    if {[regexp ".gz$" $vcfF]} {
		set f [open "|gzip -cd $vcfF"]
	    } else {
		set f [open "$vcfF"]
	    }
	    while {! [eof $f]} {
		set L [gets $f]
		if {[string range $L 0 5] eq "#CHROM"} {
		    set L_samplesOkFROMvcfF {}
		    foreach sample $g_AnnotSV(candidateSnvIndelSamples) {
			set i_sample($sample) [lsearch -exact [split $L "\t"] $sample]
			if {$i_sample($sample) ne -1} {lappend L_samplesOkFROMvcfF $sample}
		    }
		    continue
		}
		if {[string index $L 0] eq "#" || $L eq ""} {continue}
		# set the variables
		set Ls [split $L "\t"] 
		set chrom [lindex $Ls 0]
		regsub -nocase "chr" $chrom "" chrom
		set pos [lindex $Ls 1]
		set filter [lindex $Ls 6]
		set formatData [lindex $Ls 8]

		# Consider only the SNV/indel (not the SV in the VCF file)
		##########################################################
		# Example of SV: 
		# - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"
		# - Type2: "<INS>", "<DEL>", ...
		# - Type3: complex rearrangements with breakends: "G]17:1584563]"
		set ref [lindex $Ls 3]
		set alt [lindex $Ls 4]
		if {[regexp "<|\\\[|\\\]" $alt]} {continue}
		set variantLength [expr {[string length $ref]-[string length $alt]}]
		if {[expr {abs($variantLength)}]>$g_AnnotSV(SVminSize)} {continue}; # it is an SV

		# keep only variant with FILTER == PASS
		if {$g_AnnotSV(snvIndelPASS) && $filter ne "PASS"} {continue}
		
		foreach sample $L_samplesOkFROMvcfF {
		    set sampleData [lindex $Ls $i_sample($sample)]
		    # set the GT
		    set j_GT [lsearch -exact [split $formatData ":"] "GT"]
		    if {$j_GT eq -1} {continue}
		    set GT ""
		    set GTsample [lindex [split $sampleData ":"] $j_GT]
		    set GTsample [split $GTsample "/|\\|"]
		    if {[lindex $GTsample 0] ne [lindex $GTsample 1] && [lindex $GTsample 0] ne ""} {
			lappend L_htzPos($chrom,$sample) "$pos"
		    }
		}
	    }
	    close $f
	}
    }

    set textToReturn {}
    
    if {$GENEchrom eq "FULL"} {
	# No annotation for a full line (an SV can overlap several genes)
	#################################################################
	set textToReturn "$L_htzPos(FULL)"
    } else {
	# Annotation of a split line
	############################
	
	if {![info exists L_samplesFromTheHeader]} {
	    # Listing of all the samples in the VCF input files
	    regsub ".*FORMAT\t" $headerOutput "" samples
	    regsub "\tAnnotation_mode.*" $samples "" samples
	    set L_samplesFromTheHeader [split $samples "\t"]
	}
	
	# Is the SV detected in each sample from $g_AnnotSV(candidateSnvIndelSamples)?
	# Answer given in the GT of the VCF input file (but not from a bedfile)
	if {[info exists g_AnnotSV(formatColNumber)]} { ;# <=> SVinputFile is a VCF
	    set formatData [lindex $Line $g_AnnotSV(formatColNumber)]
	    set i_gt [lsearch [split $formatData ":"] "GT"]
	    foreach sample $L_samplesFromTheHeader valueSample [lrange $Line [expr {$g_AnnotSV(formatColNumber)+1}] end-13] {
		set gt [lindex [split $valueSample ":"] $i_gt]
		if {[regexp "\[1-9\]" $gt]} {
		set SV($sample) detected
		}
	    }
	}
	
	foreach sample $g_AnnotSV(candidateSnvIndelSamples) {
	    # If the SV is not detected in this sample:
	    if {![info exists SV($sample)]} {lappend textToReturn ""; continue}
	    # if $GENEchrom not present in the VCF file:
	    if {![info exists L_htzPos($GENEchrom,$sample)]} {lappend textToReturn ""; continue}
	    # List for htz variants
	    set i_first [DichotomySearch $GENEstart $L_htzPos($GENEchrom,$sample)]
	    set i_last [DichotomySearch $GENEend $L_htzPos($GENEchrom,$sample)]
	    if {$GENEstart ne [lindex $L_htzPos($GENEchrom,$sample) $i_first]} {set i_first [expr {$i_first+1}]}
	    set sampleText {}
	    foreach pos [lsort -unique [lrange $L_htzPos($GENEchrom,$sample) $i_first $i_last]] {
		lappend sampleText "${GENEchrom}_$pos"
	    }	    
	    lappend textToReturn [join $sampleText ";"]
	}
    }

    return "[join $textToReturn "\t"]" 

}
