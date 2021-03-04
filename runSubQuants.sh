#!/bin/bash
# kallisto script

for file in ./*_1.sub.fastq
do
    fname=${file##*/}
    base=${fname%_1*}
    time kallisto quant -i index.idx -o "kallisto_output/${base}.1" -b 30 -t 2 "${base}_1.sub.fastq"  "${base}_2.sub.fastq"
done