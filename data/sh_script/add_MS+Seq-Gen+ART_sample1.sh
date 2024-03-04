#!/bin/bash

for ((i=0; i<=999; i++))
do
    modified_file="./sam-ms/trisomy-ms/${i}_sample_modified.sam"
    original_file="./sam-ms/sample_${i}.sam"

    cat "$modified_file" >> "$original_file"
    echo "Appended $modified_file to $original_file"
done