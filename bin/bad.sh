#!/usr/bin/env bash
FILE="hello world.txt"
if [ -f $FILE ]; then
  echo found
fi