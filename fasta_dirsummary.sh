#!/bin/bash

# Script bash fasta_dirsummary.sh <barcode_dir>
# Autor: Luciano Kalabric
# Conta os arquivos .fastq do diretório barcode_dir

# Conta o total de arquivos .fasta do diretório
echo "Número de arquivos .fasta:"
find -H $1 -name *.fasta | wc -l

# Conta o total de reads do diretório
echo "Total reads"
grep -o ">" $(find -H  $1 -name *.fasta) | wc -l

# Conta o total de reads por BC
echo "Total reads por BC:"
for i in $(find $1/*.fasta -type f -exec basename {} .fasta \; | cut -d_ -f1 | sort); do
	grep -o ">" $1/${i}.fasta | wc -l
done
