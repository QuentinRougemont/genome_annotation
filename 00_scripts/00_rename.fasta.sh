#!/bin/bash

#very basic code that will be improved later 
source config/config

genome=$1 	#name of the genome
new_name=$2 	#name for the scaffold ids

sed "s/^>/>""$new_name""_/g" "$genome" > "$new_name".fa

