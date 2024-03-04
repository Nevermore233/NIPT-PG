#!/bin/bash

for ((i=0; i<=999; i++))
do
    input_file="./sam/trisomy/${i}_sample.sam"
    output_file="./sam/trisomy/${i}_sample_modified.sam"
    
    sed '1,3d' "$input_file" > "$output_file"
    
    echo "Processed $input_file, saved as $output_file"
done