#!/bin/bash -x

set -eo pipefail


# Run in background AnnotSV + creation of a flag ("running.flag") during the process
./commandAnnotSV.sh &

# RAM tracking
rm -f ramAnnotSV-*-ram.txt
ramFile="ramAnnotSV-`hostname`-ram.txt"

sleep 1
pid=`ps -ef | grep AnnotSV | grep -v "grep" | grep "bin/AnnotSV" | awk '{print $2}'`
while [ -e running.flag ]
do
	# "rss" = physical memory (Ko) allocated at the moment
	# The last command fails, the "echo" command will run anyway thanks to the "||"
        ps -o rss= -p $pid >> $ramFile ||  echo "End of the while"
        sleep 0.01
done

max=`sort -n $ramFile  | tail -1`
echo "rss = $max Ko"

#rm -f $ramFile


echo "ok - Finished"

