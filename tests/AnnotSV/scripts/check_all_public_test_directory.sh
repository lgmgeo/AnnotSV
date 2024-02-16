#!/bin/bash

# Timing: test_01  --> test_34 
#         13:18:21     14:30:55
# ~1h15

for f in test_*/command_Tcl.sh
do
        d=`dirname $f`
        cd $d
        echo "######## $d "
	echo "###############################################################################"
        date
        command_Tcl.sh &> command_Tcl.log
        grep -i error command_Tcl.log
        grep "ok - Finished" command_Tcl.log
        tail -2 command_Tcl.log
        echo ""
        cd ..
done


echo "ok - Finished all checkings"

