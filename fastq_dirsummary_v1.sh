#!/bin/bash

# Script bash fastq_dirsummary.sh <barcode_dir>
# Autor: Luciano Kalabric
# Conta os arquivos .fastq do diretório barcode_dir

ls $1/*.fastq | wc -l
# Conta o total de sequencias do diretório barcode_dir
grep -c "runid" $1/*.fastq | cut -d : -f 2 | awk '{s+=$1} END {printf "%.0f\n",s}'

