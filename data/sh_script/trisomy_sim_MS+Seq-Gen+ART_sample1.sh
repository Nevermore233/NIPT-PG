#!/bin/bash
#illumina test examples
art=../art_illumina
for ((i=0; i<=333; i++))
do
    output_file="${i}_sample"
 
    $art -ss HS25 -sam -i ./ref_MS+Seq-Gen+ART_sample1_chr13.fa -l 35 -p -f 0.2  -m 100 -s 5 -M -ir 0.005 -ir2 0.005 -dr 0.005 -dr2 0.005 -na -o "./sam-ms/trisomy-ms/${output_file}"
done

for ((i=334; i<=666; i++))
do
    output_file="${i}_sample"
 
    $art -ss HS25 -sam -i ./ref_MS+Seq-Gen+ART_sample1_chr18.fa -l 35 -p -f 0.2  -m 100 -s 5 -M -ir 0.005 -ir2 0.005 -dr 0.005 -dr2 0.005 -na -o "./sam-ms/trisomy-ms/${output_file}"
done

for ((i=667; i<=999; i++))
do
    output_file="${i}_sample"
 
    $art -ss HS25 -sam -i ./ref_MS+Seq-Gen+ART_sample1_chr21.fa -l 35 -p -f 0.2  -m 100 -s 5 -M -ir 0.005 -ir2 0.005 -dr 0.005 -dr2 0.005 -na -o "./sam-ms/trisomy-ms/${output_file}"
done