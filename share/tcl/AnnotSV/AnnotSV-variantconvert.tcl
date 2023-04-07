############################################################################################################
# AnnotSV 3.3.3                                                                                            #
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
# - Check if the "pip install -e ." command was already run
proc checkVariantconvertConfigfile {} {

    global g_AnnotSV

    # Use of variantconvert to create a vcf output
    if {$g_AnnotSV(vcf)} {
	
	# - Check if the reference fasta file in the variantconvert bed configfiles has the good path
        #############################################################################################

	## SVinputfile is a BED
	## (no need to have a reference fasta file from a VCF SVinputfile)
	if {[regexp "\\.bed$" $g_AnnotSV(SVinputFile)]} {

            set configfile "$g_AnnotSV(variantconvertDir)/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_bed.json"
            set localconfigfile "$g_AnnotSV(variantconvertDir)/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_bed.local.json"

	    set distributedPathLine "\"path\": \".*\","
	    set newPathLine         "\"path\": \"$g_AnnotSV(installDir)/share/AnnotSV/Annotations_Human/BreakpointsAnnotations/GCcontent/GRCh37/GRCh37_chromFa.fasta\","
	    set distributedRefLine  "\"##reference=file:.*\""
	    set newRefLine          "\"##reference=file:$g_AnnotSV(installDir)/share/AnnotSV/Annotations_Human/BreakpointsAnnotations/GCcontent/GRCh37/GRCh37_chromFa.fasta\""

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
	
        # - Check if the "pip install -e ." command was already run
        ###########################################################
        if {![file exists $g_AnnotSV(variantconvertDir)/pipinstall.flag]} {
            set currentDir [pwd]
            cd $g_AnnotSV(variantconvertDir)
            if {[catch {exec pip3 install -e .}]} {
                if {[catch {exec pip install -e .} Message]} {
                    WriteTextInFile "$Message" $g_AnnotSV(variantconvertDir)/pipinstall.flag
                } else {
                    WriteTextInFile "Done" $g_AnnotSV(variantconvertDir)/pipinstall.flag
                }
            } else {
                WriteTextInFile "Done" $g_AnnotSV(variantconvertDir)/pipinstall.flag
            }
            cd $currentDir
        }
    }
}




