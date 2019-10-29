#!/bin/bash
ss -tln |
  awk 'NR > 1{gsub(/.*:/,"",$4); print $4}' |
  sort -un |
  awk -v n=50000 '$0 < n {next}; $0 == n {n++; next}; {exit}; END {print n}'

