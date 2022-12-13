#!/bin/bash

# com: Compare two sorted files line by line. Output the lines that are common, plus the lines that are unique.
#      comm [options]... File1 File2
#       -2 Suppress lines unique to file2
#       -3 Suppress lines that appear in both files
#
# seq: print a sorted sequence of numbers from 5000 to 50000
#
# ss: used to dump socket statistics.
#       -H, --no-header
#       -n, --numeric
#              Do not try to resolve service names. Show exact bandwidth
#              values, instead of human-readable.
#       -t, --tcp
#              Display TCP sockets.
#     e.g. "ss -Htan" >>> "LISTEN    0      5                       0.0.0.0:63864                  0.0.0.0:*"

comm -23 <(seq 5000 50000 | sort) <(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) | head -n 1
