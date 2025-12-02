#!/bin/sh

# USAGE:
# bash ./check_catalogue.sh *_merge.pksc.*
catalogue=$1

echo "${catalogue}"

# Link and compile
gfortran check_catalogue.f90 -o check_catalogue

# Execute
./check_catalogue << EOF
'${catalogue}'
EOF

# Clean up
rm -r -f check_catalogue.mod check_catalogue.o
