#!/bin/bash

comm -23 <(seq 5000 50000 | sort) <(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) | head -n 1
