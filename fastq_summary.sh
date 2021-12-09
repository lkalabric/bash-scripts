#!/bin/bash

# Script bash fastq_summary.sh <fastq_dir> ou <fastq file>
# Autor: Luciano Kalabric Conta
# Conta os arquivos .fastq de um diretório e as readas dos arquivos

if [[ $# -eq 0 ]]; then
	echo "Falta o nome do diretório ou arquivo .fastq!"
	exit 0
fi
if [[ -d $1 ]]; then
	# Conta os arquivos do diretório fastq_dir
	ls $1/*.fastq | wc -l
	# Conta o total de reads de todos os arquivos contidos no diretório fastq_dir
	grep -c "runid" $1/*.fastq | cut -d : -f 2 | awk '{s+=$1} END {printf "%.0f\n",s}'
elif [[ -f $1 ]]; then
	# Conta o total de reads do arquivo .fastq
	grep -c "runid" $1 | cut -d : -f 2 | awk '{s+=$1} END {printf "%.0f\n",s}'
else
	echo "Diretório ou arquivo inexistente."
fi
