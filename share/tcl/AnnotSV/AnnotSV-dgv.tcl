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

## - Check if the DGV files have been downloaded:
#    - DGV.GS.*.gff3
#    - *_supportingvariants_*.txt
#
## - Check and create if necessary the following files:
#    - 'date'_DGV.GS*_annotations.sorted.bed
#    - 'date'_DGV_samplesInStudies.tsv
proc checkDGVfiles {} {

    global g_AnnotSV

    ## Check if the 2 DGV files have been downloaded then formatted
    ##############################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set DGVfile1Downloaded [glob -nocomplain "$extannDir/SVincludedInFt/DGV/$g_AnnotSV(genomeBuild)/DGV.GS.*.gff3"]
    set DGVfile2Downloaded [glob -nocomplain "$extannDir/SVincludedInFt/DGV/$g_AnnotSV(genomeBuild)/*_supportingvariants_*.txt"]
    set DGVfile1Formatted [glob -nocomplain "$extannDir/SVincludedInFt/DGV/$g_AnnotSV(genomeBuild)/*_DGV.GS*_annotations.sorted.bed"]
    set DGVfile2Formatted [glob -nocomplain "$extannDir/SVincludedInFt/DGV/$g_AnnotSV(genomeBuild)/*_DGV_samplesInStudies.tsv"]

    if {($DGVfile1Downloaded ne "" && $DGVfile2Downloaded ne "") || ($DGVfile1Formatted ne "" && $DGVfile2Formatted ne "")} {
	# DGV annotation
	set g_AnnotSV(dgvAnn) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "DGV_GAIN_IDs DGV_GAIN_n_samples_with_SV DGV_GAIN_n_samples_tested DGV_GAIN_Frequency DGV_LOSS_IDs DGV_LOSS_n_samples_with_SV DGV_LOSS_n_samples_tested DGV_LOSS_Frequency" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(dgvAnn) 0; return}
    } else {
	# No DGV annotation
	set g_AnnotSV(dgvAnn) 0
	return
    }

    if {[llength $DGVfile1Formatted]>1} {
	puts "Several DGV files exist:"
	puts "$DGVfile1Formatted"
	puts "Keep only one: [lindex $DGVfile1Formatted end]\n"
	foreach dgv [lrange $DGVfile1Formatted 0 end-1] {
	    file rename -force $dgv $dgv.notused
	}
	return
    }
    if {[llength $DGVfile2Formatted]>1} {
	puts "Several DGV files exist:"
	puts "$DGVfile2Formatted"
	puts "Keep only one: [lindex $DGVfile2Formatted end]\n"
	foreach dgv [lrange $DGVfile2Formatted 0 end-1] {
	    file rename -force $dgv $dgv.notused
	}
	return
    }

    if {$DGVfile1Formatted eq "" || $DGVfile2Formatted eq ""} {
	# The downloaded file exist but not the 2 formatted.
	#
	## Create:
	##   - 'date'_DGV.GS_annotations.sorted.bed file     ; # Header: chr start end variantsubtype supportingsamples studiessamples
	##   - 'date'_DGV_samplesInStudies.tsv               ; # Syntax: Study num_unique_samples_tested {samples}
	#######################################################################################################################

	## Create : 'date'_DGV.GS_annotations.bed
	set DGVfile1Formatted "$extannDir/SVincludedInFt/DGV/$g_AnnotSV(genomeBuild)/[clock format [clock seconds] -format "%Y%m%d"]_DGV.GS_annotations.sorted.bed"

	puts "...DGV configuration \[part 1\] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	puts "\t...creation of $DGVfile1Formatted"
	puts "\t   (done only once during the first DGV annotation)"

	set f [open "$DGVfile1Downloaded"]
	while {! [eof $f]} {
	    set L [gets $f]
	    set Ls [split $L "\t"]
	    set infos [split [lindex $Ls 8] ";"]
	    set i_ID             [lsearch -regexp $infos "^ID"]; if {$i_ID == -1} {puts "Bad syntax into $DGVfile1Downloaded.\nID field not found - Exit with error"; exit 2}
	    set i_variantsubtype [lsearch -regexp $infos "^variant_sub_type"]; if {$i_variantsubtype == -1} {puts "Bad syntax into $DGVfile1Downloaded.\nvariant_sub_type field not found - Exit with error"; exit 2}
	    set i_start          [lsearch -regexp $infos "^outer_start"]; if {$i_start == -1} {puts "Bad syntax into $DGVfile1Downloaded.\nouter_start field not found - Exit with error"; exit 2}
	    set i_end            [lsearch -regexp $infos "^outer_end"]; if {$i_end == -1} {puts "Bad syntax into $DGVfile1Downloaded.\nouter_end field not found - Exit with error"; exit 2}
	    set i_variants       [lsearch -regexp $infos "^variants"]; if {$i_variants == -1} {puts "Bad syntax into $DGVfile1Downloaded.\nvariants field not found - Exit with error"; exit 2}
	    set i_studies        [lsearch -regexp $infos "^Studies"]; if {$i_studies == -1} {puts "Bad syntax into $DGVfile1Downloaded.\nStudies field not found - Exit with error"; exit 2}
	    set i_samples        [lsearch -regexp $infos "^samples"]; if {$i_samples == -1} {puts "Bad syntax into $DGVfile1Downloaded.\nsamples field not found - Exit with error"; exit 2}
	    break
	}
	close $f
	set TexteToWrite ""
	set i -1
	foreach L [LinesFromFile $DGVfile1Downloaded] {
	    incr i
	    if {[eval expr {$i%3}]} {continue}
	    set Ls [split $L "\t"]

	    # 3 lines for each ID::
	    # chr1    SV     copy_number_variation_region    812998  812998  .       .       .       ID=gssvG1;Name=gssvG1;variant_type=SV;variant_sub_type=Gain;outer_start=812998;inner_start=837847;inner_end=1477469;outer_end=1708649;inner_rank=10;num_variants=3;variants=nssv24621,nssv1423341,nssv1440790;num_studies=2;Studies=Perry2008,Park2010;num_platforms=2;Platforms=AgilentCustom_015685+015686+244K,Agilent24M;number_of_algorithms=1;algorithms=ADM2;num_samples=3;samples=NA18968,NA18969,NA19221;Frequency=5.45%;PopulationSummary=African 1:Asian 2:European 0:Mexican 0:MiddleEast 0:NativeAmerican 0:NorthAmerican 0:Oceania 0:SouthAmerican 0:Turkish 0:Admixed 0:Unknown 0;num_unique_samples_tested=55
	    # chr1    SV     copy_number_variation_region    837847  1477469 .       .       .       ID=gssvG1;Name=gssvG1;variant_type=SV;variant_sub_type=Gain;outer_start=812998;inner_start=837847;inner_end=1477469;outer_end=1708649;inner_rank=10;num_variants=3;variants=nssv24621,nssv1423341,nssv1440790;num_studies=2;Studies=Perry2008,Park2010;num_platforms=2;Platforms=AgilentCustom_015685+015686+244K,Agilent24M;number_of_algorithms=1;algorithms=ADM2;num_samples=3;samples=NA18968,NA18969,NA19221;Frequency=5.45%;PopulationSummary=African 1:Asian 2:European 0:Mexican 0:MiddleEast 0:NativeAmerican 0:NorthAmerican 0:Oceania 0:SouthAmerican 0:Turkish 0:Admixed 0:Unknown 0;num_unique_samples_tested=55
	    # chr1    SV     copy_number_variation_region    1708649 1708649 .       .       .       ID=gssvG1;Name=gssvG1;variant_type=SV;variant_sub_type=Gain;outer_start=812998;inner_start=837847;inner_end=1477469;outer_end=1708649;inner_rank=10;num_variants=3;variants=nssv24621,nssv1423341,nssv1440790;num_studies=2;Studies=Perry2008,Park2010;num_platforms=2;Platforms=AgilentCustom_015685+015686+244K,Agilent24M;number_of_algorithms=1;algorithms=ADM2;num_samples=3;samples=NA18968,NA18969,NA19221;Frequency=5.45%;PopulationSummary=African 1:Asian 2:European 0:Mexican 0:MiddleEast 0:NativeAmerican 0:NorthAmerican 0:Oceania 0:SouthAmerican 0:Turkish 0:Admixed 0:Unknown 0;num_unique_samples_tested=55

	    regsub "chr" [lindex $Ls 0] "" chr

	    set infos [split [lindex $Ls 8] ";"]

	    set ID [lindex $infos $i_ID];                         regsub ".*=" $ID "" ID
	    set variantsubtype [lindex $infos $i_variantsubtype]; regsub ".*=" $variantsubtype "" variantsubtype
	    set start [lindex $infos $i_start];                   regsub ".*=" $start "" start
	    set end [lindex $infos $i_end];                       regsub ".*=" $end "" end
	    set variants [lindex $infos $i_variants];             regsub ".*=" $variants "" variants
	    set studies [lindex $infos $i_studies];               regsub ".*=" $studies "" studies
	    set samples [lindex $infos $i_samples];               regsub ".*=" $samples "" samples

	    lappend L_studFromGS {*}[split $studies ","]
	    lappend TexteToWrite "$chr\t$start\t$end\t$ID\t$variantsubtype\t$studies\t$variants\t$samples"
	}
	WriteTextInFile [join $TexteToWrite "\n"] $DGVfile1Formatted.tmp
	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
	if {[catch {eval exec sort -k1,1 -k2,2n $DGVfile1Formatted.tmp > $DGVfile1Formatted} Message]} {
	    puts "-- checkDGVfiles --"
	    puts "sort -k1,1 -k2,2n $DGVfile1Formatted.tmp > $DGVfile1Formatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $DGVfile1Formatted.tmp

	set L_studFromGS [lsort -unique $L_studFromGS]


	## Create : 'date'_DGV_samplesInStudies.tsv
	set DGVfile2Formatted "$extannDir/SVincludedInFt/DGV/$g_AnnotSV(genomeBuild)/[clock format [clock seconds] -format "%Y%m%d"]_DGV_samplesInStudies.tsv"

	puts "...DGV configuration \[part 2\] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	puts "\t...creation of $DGVfile2Formatted"
	puts "\t   (done only once during the first DGV annotation)\n"

	# Parsing of $DGVfile2Downloaded (=file with all supporting variants (not merged))
	foreach L [LinesFromFile $DGVfile2Downloaded] {
	    set Ls [split $L "\t"]

	    # Header + example line:
	    ########################
	    # variantaccession        chr     start   end     varianttype     variantsubtype  reference       pubmedid        method  platform        mergedvariants  supportingvariants      mergedorsample  frequency       samplesize  observedgains   observedlosses  cohortdescription       genes   samples
	    # nssv1750706     1       10001   19818   SV     duplication     Sudmant_et_al_2013      23825009        Oligo aCGH,Sequencing           nsv945697       ""      S               97      1       0               DDX11L1,MIR6859-1,MIR6859-2,WASH7P  HGDP00521
	    if {[regexp "^variantaccession" $L]} {
		set i_reference [lsearch -exact $Ls "reference"];   if {$i_reference == -1} {puts "$DGVfile2Downloaded"; puts "Bad header line syntax. reference column not found - Exit with error"; exit 2}
		set i_samples  [lsearch -exact $Ls "samples"];    if {$i_samples == -1} {puts "$DGVfile2Downloaded"; puts "Bad header line syntax. samples column not found - Exit with error"; exit 2}
		set i_samplesize  [lsearch -exact $Ls "samplesize"];    if {$i_samplesize == -1} {puts "$DGVfile2Downloaded"; puts "Bad header line syntax. samplesize column not found - Exit with error"; exit 2}
		continue
	    }

	    # NOTE: only 1 reference by variant (only variants from a same study have been merged)

	    # WARNING: Some merged variants (esv or nsv) are associated with several samples:
	    #   -> esv2740664: "SSM008,SSM017,SSM053".
	    # These sv correspond to ssv present in this file, merged from the same study.
	    # => No needs to parse
	    #
	    # WARNING: Some studies have only sv (and no ssv!!) (ex: Kidd_et_al_2010b, Boomsma_et_al_2014, Altshuler_et_al_2010, ...)
	    #          => We have to consider the sv (and not only the ssv)
	    set reference  [lindex $Ls $i_reference]
	    regsub -all "_et_al_|-|_|Consortium" $reference "" reference ; # To have the same format as in DGV.GS
	    set samples    [lindex $Ls $i_samples]
	    set samplesize [lindex $Ls $i_samplesize]
	    # WARNING: Some variants are not associated with a sample:
	    #   -> GRCh37_hg19_supportingvariants_2016-05-15.txt: 496 427 variants not associated / 6 668 716 variants
	    set size($reference) $samplesize
	    lappend lSamples($reference) {*}[split $samples ","]
	}

	## Correction of inconsistence in DGV march 2016 (given by J. MacDonald, march 2017)
	set size(Durbin2010) $size(1000GenomesPilotProject) ; # Durbin2010 = 1000GenomesPilotProject (only the tandem duplications)
	set lSamples(Durbin2010) $lSamples(1000GenomesPilotProject)
	set size(Arlt2010) $size(Arlt2011) ; # Arlt2010 = Arlt2011
	set lSamples(Arlt2010) $lSamples(Arlt2011)
	set size(Coe2014) 11256 ; # errors in the submission: report of the numbers for cases and controls even though the case SV data wasn’t included in DGV
	set size(Cooper2011) 8329 ; # errors in the submission: report of the numbers for cases and controls even though the case SV data wasn’t included in DGV


	foreach study [lsort [array names lSamples]] {
	    set lSamples($study) [lsort -unique $lSamples($study)]
	    WriteTextInFile "$study\t$size($study)\t$lSamples($study)" $DGVfile2Formatted
	}

	# Is there any study from DGV.GS that are not in the DGV primary data?
	# DGV march 2016:
	#    2 studies from GS not in DGV primary data:
	#    Arlt2010 and Durbin2010 !!!
	foreach study $L_studFromGS {
	    if {![info exists lSamples($study)]} {
		puts $study
	    }
	}
	puts "\n"

	file delete -force $DGVfile1Downloaded
	file delete -force $DGVfile2Downloaded
    }

    if {$DGVfile1Formatted eq "" || $DGVfile2Formatted eq ""} {
	puts "No DGV annotation available"
	set g_AnnotSV(dgvAnn) 0
    }
}



proc DGVannotation {SVchrom SVstart SVend L_i} {

    global g_AnnotSV
    global dgvText


    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set DGVfile1Formatted [glob -nocomplain "$extannDir/SVincludedInFt/DGV/$g_AnnotSV(genomeBuild)/*_DGV.GS*_annotations.sorted.bed"]
    set DGVfile2Formatted [glob -nocomplain "$extannDir/SVincludedInFt/DGV/$g_AnnotSV(genomeBuild)/*_DGV_samplesInStudies.tsv"]


    if {![info exists dgvText(DONE)]} {

	# headerOutput "\tDGV_GAIN_IDs\tDGV_GAIN_n_samples_with_SV\tDGV_GAIN_n_samples_tested\tDGV_GAIN_Frequency"
	# headerOutput "\tDGV_LOSS_IDs\tDGV_LOSS_n_samples_with_SV\tDGV_LOSS_n_samples_tested\tDGV_LOSS_Frequency"
	set L_dgvText(Empty) "{} {0} {0} {-1} {} {0} {0} {-1}"
	foreach i $L_i {
	    lappend dgvText(Empty) "[lindex $L_dgvText(Empty) $i]"
	}
	set dgvText(Empty) [join $dgvText(Empty) "\t"]

	# Loading studies information
	foreach L [LinesFromFile $DGVfile2Formatted] {
	    set Ls [split $L "\t"]
	    set study [lindex $Ls 0]
	    set number($study) [lindex $Ls 1]
	    set samples($study) [lindex $Ls 2]
	}

	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.dgv" tmpFile
	set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	file delete -force $tmpFile
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted  -a $g_AnnotSV(fullAndSplitBedFile) -b $DGVfile1Formatted -wa -wb > $tmpFile} Message]} {
	    puts "-- DGVannotation --"
	    puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $DGVfile1Formatted -wa -wb > $tmpFile"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}

	# Parse
	set f [open $tmpFile]
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    set Ls [split $L "\t"]
	    set SVtoAnn_chrom [lindex $Ls 0]
	    set SVtoAnn_start [lindex $Ls 1]
	    set SVtoAnn_end   [lindex $Ls 2]

	    set SVdgv_chrom    [lindex $Ls end-7]
	    set SVdgv_start    [lindex $Ls end-6]
	    set SVdgv_end      [lindex $Ls end-5]
	    set SVdgv_ID       [lindex $Ls end-4]
	    set SVdgv_type     [lindex $Ls end-3]
	    set SVdgv_studies  [lindex $Ls end-2]
	    set SVdgv_variants [lindex $Ls end-1]
	    set SVdgv_samples  [lindex $Ls end]

	    # Select:
	    # - DGV SV that share > XX % length with the SV to annotate
	    # - insertion of DGV SV inside the SV to annotate
	    set SVdgv_length [expr {$SVdgv_end - $SVdgv_start}]
	    set SVtoAnn_length [expr {$SVtoAnn_end - $SVtoAnn_start}]
	    # The dgv SV is an insertion or a breakpoint
	    if {$SVdgv_length<1} {
		set SVdgv_length 1
	    }
	    # The SV to annotate is an insertion or a breakpoint
	    if {$SVtoAnn_length<1} {
		set SVtoAnn_length 1
	    }

	    if {$SVtoAnn_start < $SVdgv_start} {
		set overlap_start $SVdgv_start
	    } else {
		set overlap_start $SVtoAnn_start
	    }
	    if {$SVtoAnn_end < $SVdgv_end} {
		set overlap_end $SVtoAnn_end
	    } else {
		set overlap_end $SVdgv_end
	    }
	    set overlap_length [expr {$overlap_end - $overlap_start}]

	    # Keeping only DGV respecting the overlaps (reciprocal or not reciprocal)
	    if {[expr {$overlap_length*100.0/$SVtoAnn_length}] < $g_AnnotSV(overlap)} {continue}
	    if {$g_AnnotSV(reciprocal) eq "yes"} {
		if {[expr {$overlap_length*100.0/$SVdgv_length}] < $g_AnnotSV(overlap)} {continue}
	    }

	    if {![info exists remember($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end)]} {
		set DGVgainID($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {}
		set DGVgainVariants($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {}
		set DGVgainSamples($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {}
		set DGVgainStudies($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {}
		set DGVlossID($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {}
		set DGVlossVariants($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {}
		set DGVlossSamples($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {}
		set DGVlossStudies($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {}

	    }
	    set remember($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) 1

	    if {$SVdgv_type eq "Gain"} {
		lappend DGVgainID($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) $SVdgv_ID
		lappend DGVgainVariants($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $SVdgv_variants ","]
		lappend DGVgainSamples($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $SVdgv_samples ","]
		lappend DGVgainStudies($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $SVdgv_studies ","]

	    } elseif {$SVdgv_type eq "Loss"} {
		lappend DGVlossID($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) $SVdgv_ID
		lappend DGVlossVariants($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $SVdgv_variants ","]
		lappend DGVlossSamples($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $SVdgv_samples ","]
		lappend DGVlossStudies($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $SVdgv_studies	","]
	    }
	}


	# Loading DGV final annotation for each SV
	# headerOutput "\tDGV_GAIN_IDs\tDGV_GAIN_n_samples_with_SV\tDGV_GAIN_n_samples_tested\tDGV_GAIN_Frequency"
	#              "\tDGV_LOSS_IDs\tDGV_LOSS_n_samples_with_SV\tDGV_LOSS_n_samples_tested\tDGV_LOSS_Frequency"
	foreach SVtoAnn [array names remember] {
	    set L_dgvText($SVtoAnn) ""
	    # DGV_GAIN_IDs
	    lappend L_dgvText($SVtoAnn) [join $DGVgainID($SVtoAnn) ","]

	    # DGV_GAIN_n_samples_with_SV
	    set lengthVariantsSorted [llength [lsort -unique $DGVgainVariants($SVtoAnn)]]
	    set lengthSamples [llength $DGVgainSamples($SVtoAnn)]
	    set lengthSamplesSorted [llength [lsort -unique $DGVgainSamples($SVtoAnn)]]
	    set doublon [expr {$lengthSamples-$lengthSamplesSorted}]
	    set n_samples_with_SV [expr {$lengthVariantsSorted-$doublon}]
	    lappend L_dgvText($SVtoAnn) "$n_samples_with_SV"

	    # DGV_GAIN_n_samples_tested
	    set allStudies [lsort -unique $DGVgainStudies($SVtoAnn)]
	    if {$allStudies ne ""} {
		set n 0
		set samp ""
		foreach study $allStudies {
		    incr n $number($study)
		    lappend samp {*}$samples($study)
		}
		set doublon [expr [llength $samp]-[llength [lsort -unique $samp]]]
		set n_samples_tested [expr {$n -$doublon}]
		lappend L_dgvText($SVtoAnn) "$n_samples_tested"

		# DGV_GAIN_Frequency
		set freq [expr {$n_samples_with_SV*1.0/$n_samples_tested}]
		set freq [format "%.8f" $freq]
		if {[set g_AnnotSV(metrics)] eq "fr"} {regsub {\.} $freq "," freq}
		lappend L_dgvText($SVtoAnn) "$freq"
	    } else {
		lappend L_dgvText($SVtoAnn) {*}{0 -1}
	    }

	    # DGV_LOSS_IDs
	    lappend L_dgvText($SVtoAnn) "[join $DGVlossID($SVtoAnn) ","]"

	    # DGV_LOSS_n_samples_with_SV
	    set lengthVariantsSorted [llength [lsort -unique $DGVlossVariants($SVtoAnn)]]
	    set lengthSamples [llength $DGVlossSamples($SVtoAnn)]
	    set lengthSamplesSorted [llength [lsort -unique $DGVlossSamples($SVtoAnn)]]
	    set doublon [expr {$lengthSamples-$lengthSamplesSorted}]
	    set n_samples_with_SV [expr {$lengthVariantsSorted-$doublon}]
	    lappend L_dgvText($SVtoAnn) "$n_samples_with_SV"

	    # DGV_LOSS_n_samples_tested
	    set allStudies [lsort -unique $DGVlossStudies($SVtoAnn)]
	    if {$allStudies ne ""} {
		set n 0
		set samp ""
		foreach study $allStudies {
		    incr n $number($study)
		    lappend samp $samples($study)
		}
		set doublon [expr [llength $samp]-[llength [lsort -unique $samp]]]

		set n_samples_tested [expr {$n -$doublon}]
		lappend L_dgvText($SVtoAnn) "$n_samples_tested"

		# DGV_LOSS_Frequency
		set freq [expr {$n_samples_with_SV*1.0/$n_samples_tested}]
		set freq [format "%.8f" $freq]
		# Change metrics from "." to ","
		if {[set g_AnnotSV(metrics)] eq "fr"} {regsub -all {\.} $freq "," freq}
		lappend L_dgvText($SVtoAnn) "$freq"
	    } else {
		lappend L_dgvText($SVtoAnn) {*}{0 -1}
	    }

	    # Keep only the user requested columns (defined in the configfile)
	    set dgvText($SVtoAnn) ""
	    foreach i $L_i {
		lappend dgvText($SVtoAnn) "[lindex $L_dgvText($SVtoAnn) $i]"
	    }
	    set dgvText($SVtoAnn) [join $dgvText($SVtoAnn) "\t"]

	}
	catch {unset L_dgvText}
	set dgvText(DONE) 1
	file delete -force $tmpFile
    }

    if {[info exist dgvText($SVchrom,$SVstart,$SVend)]} {
	return $dgvText($SVchrom,$SVstart,$SVend)
    } else {
	return $dgvText(Empty)
    }
}
