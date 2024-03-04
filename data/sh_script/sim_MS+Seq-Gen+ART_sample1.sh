#!/bin/bash
#illumina test examples
art=../art_illumina
f_value="f_value_MS+Seq-Gen+ART_sample1.txt"
> $f_value

for ((i=0; i<=1999; i++))
do
    output_file="sample_${i}"
    
    # Generate a random floating point number within the range of [0.1, 0.5]
    random_float=$(awk -v min=0.1 -v max=0.5 'BEGIN{srand(); print min+rand()*(max-min)}')
    
    echo $random_float >> $f_value
     
    $art -ss HS25 -sam -i ./ref_MS+Seq-Gen+ART_sample1.fa -l 35 -p -f $random_float -m 100 -s 5 -M -ir 0.005 -ir2 0.005 -dr 0.005 -dr2 0.005 -na -o "./sam-ms/${output_file}"
done