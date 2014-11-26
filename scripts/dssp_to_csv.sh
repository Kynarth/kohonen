#!/usr/bin/env bash

#----------------------------------------------------------------------------#
# From dssp file, we get only the informations about each AA from a protein. #
# The result, in csv format, is added to a mutual csv file.                  #
#----------------------------------------------------------------------------#
file=$1

# Create result directory if not exist
if [ ! -d "results" ]; then
    mkdir results
fi

# Get protein id
id=$(basename $file .dssp)

# Get the header
line=$(grep -n "#" $file | cut -d":" -f1)

# Delete header and useless informations
sed "1,${line}d" $file | \

# Replace tab spliting by ';'
awk -v FIELDWIDTHS='5 7 3 3 1 1 1 1 1 1 1 4 4 1 4 12 11 11 11 8 6 6 6 6 7 7 7' \
    -v OFS=';' '{ $1=$1 ""; print }' | \

# Replace empty structure assign by coil > 'C'
awk 'BEGIN{OFS=FS=";"}$4=="   "{$4="  C"}{print}' | \

# Delete spaces
sed -r 's/\s+//g' | \

# Add protein id and delete 'b' AA
sed "s/^/$id;/g" | awk -F";" '$4 != "b"' >> results/mutual.csv

# Add string delimiter and append it in mutual csv file
# sed -r 's/(-?[0-9]*(([a-Z]+[0-9]*)|[-<+>!*]|-?[0-9],-?[0-9].[0-9])+);/"\1";/g' >> results/mutual.csv
