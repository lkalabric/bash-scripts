#!/bin/bash

# Script bash fastq_dirsummary.sh <barcode_dir>
# Autor: Luciano Kalabric
# Conta os arquivos .fastq do diretório barcode_dir

find -H $1 -name *.fastq | wc -l
# Conta o total de sequencias do diretório barcode_dir
grep -o "runid" $(find -H  $1 -name *.fastq) | wc -l
# grep -c "runid" $1/*.fastq | cut -d : -f 2 | awk '{s+=$1} END {printf "%.0f\n",s}'

