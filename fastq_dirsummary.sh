#!/bin/bash

# Script bash fastq_dirsummary.sh <barcode_dir>
# Autor: Luciano Kalabric
# Conta os arquivos .fastq do diretório barcode_dir

# Conta o total de arquivos .fastq do diretório
echo "Número de arquivos .fastq:"
find -H $1 -name *.fastq | wc -l

# Conta o total de reads do diretório
echo "Total reads"
grep -o "runid" $(find -H  $1 -name *.fastq) | wc -l
# grep -c "runid" $1/*.fastq | cut -d : -f 2 | awk '{s+=$1} END {printf "%.0f\n",s}'

# Conta o total de reads por BC
echo "Total reads por BC:"
for i in $(find $1/*.fastq -type f -exec basename {} .fastq \; | cut -d_ -f1 | sort); do
	grep -o "runid" $1/${i}.fastq | wc -l
done
