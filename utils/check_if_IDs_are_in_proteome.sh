#!/bin/bash

database=$1
ID_list=$2

while read line; do
	for ID in ${line[@]}; do
		grep "^>$ID" $database > /dev/null  || echo "Missing ID: $ID"
	done
done < "$ID_list"

