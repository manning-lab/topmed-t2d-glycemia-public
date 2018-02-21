#!/bin/bash

if [ ! -f duplicates.txt ]
then
    Rscript duplicates.R
fi

while [ ! -f duplicates.txt ]
do
    echo "waiting for duplicates.txt"
    sleep 1
done

if [[ -z "$1" ]]
then
    echo "Enter trait name from phenotype file"
    exit 1
else 
    echo "trait is $1"
	mkdir $1
    awk -v t=$1 -v p="TRAIT" '{gsub(p,t)}1' duplicates.R > remove_duplicates.$1.R
    Rscript remove_duplicates.$1.R $1
fi

exit


