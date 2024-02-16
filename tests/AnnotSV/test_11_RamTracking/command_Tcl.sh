#!/bin/bash -x

set -eo pipefail


# Run AnnotSV + creation of a flag ("running.flag") during the process
./commandAnnotSV.sh &

# RAM tracking
rm -f ramAnnotSV-`hostname`-ram.txt
sleep 1
pid=`ps -ef | grep AnnotSV | grep -v "grep" | grep "bin/AnnotSV" | awk '{print $2}'`
while [ -e running.flag ]
do
	# "rss" = physical memory (Ko) allocated at the moment
	# The last command fails, the "echo" command will run anyway thanks to the "||"
        ps -o rss= -p $pid >> ramAnnotSV-`hostname`-ram.txt ||  echo "End of the while"
        sleep 0.01
done

ramFile="ramAnnotSV-`hostname`-ram.txt"
max=`sort -n $ramFile  | tail -1`
echo "rss = $max Ko"

echo "ok - Finished"

