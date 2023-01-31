############################################################################################################
# AnnotSV 3.2.3                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2023 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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


# - Check if the reference fasta file in the variantconvert bed configfiles has the good path
proc checkVariantconvertConfigfile {} {

    global g_AnnotSV
    global env 

    # Use of variantconvert to create a vcf output
    if {$g_AnnotSV(vcf)} {
	
        ## SVinputfile is a BED
	## (no need to have a reference fasta file from a VCF SVinputfile)
	if {[regexp "\\.bed$" $g_AnnotSV(SVinputFile)]} {

            set g_AnnotSV(pythonDir) "$g_AnnotSV(installDir)/share/python3"
	    set variantconvertDIR "$g_AnnotSV(pythonDir)/variantconvert"

	    if {$g_AnnotSV(genomeBuild) == "GRCh37"} {
                set configfile "$variantconvertDIR/configs/config_annotsv3_from_bed_GRCh37.json"
                set localconfigfile "$variantconvertDIR/configs/config_annotsv3_from_bed_GRCh37.local.json"

		set distributedPathLine "\"path\": \"human_g1k_v37.fasta\""
		set newPathLine         "\"path\": \"$env(ANNOTSV)/share/AnnotSV/Annotations_Human/BreakpointsAnnotations/GCcontent/GRCh37/GRCh37_chromFa.fasta\""
		set distributedRefLine  "##reference=file:human_g1k_v37.fasta"
		set newRefLine          "##reference=file:$env(ANNOTSV)/share/AnnotSV/Annotations_Human/BreakpointsAnnotations/GCcontent/GRCh37/GRCh37_chromFa.fasta"

            } elseif {$g_AnnotSV(genomeBuild) == "GRCh38"} {
                set configfile "$variantconvertDIR/configs/config_annotsv3_from_bed_GRCh38.json"
                set localconfigfile "$variantconvertDIR/configs/config_annotsv3_from_bed_GRCh38.local.json"
		
		set distributedPathLine "\"path\": \"Homo_sapiens_assembly38.fasta\""
                set newPathLine "\"path\": \"$env(ANNOTSV)/share/AnnotSV/Annotations_Human/BreakpointsAnnotations/GCcontent/GRCh38/GRCh38_chromFa.fasta\""
		set distributedRefLine  "##reference=file:Homo_sapiens_assembly38.fasta"
		set newRefLine          "##reference=file:$env(ANNOTSV)/share/AnnotSV/Annotations_Human/BreakpointsAnnotations/GCcontent/GRCh38/GRCh38_chromFa.fasta"
            }
	    
	    if {![file exists $localconfigfile]} {
		set L_Lines {}
		foreach L [LinesFromFile $configfile] {
		    if {[regexp "$distributedPathLine" $L]} {
			regsub "$distributedPathLine" $L "$newPathLine" L
		    }
		    if {[regexp "$distributedRefLine" $L]} {
			regsub "$distributedRefLine" $L "$newRefLine" L
		    }
		    lappend L_Lines $L
		}
		ReplaceTextInFile [join $L_Lines "\n"] $localconfigfile
	    }
	}
    }
}

