#!/bin/bash

# Script bash fast5_dirsummary.sh <fast5_dir>
# Autor: Luciano Kalabric
# Conta os arquivos .fastq do diretório fast5_dir

# Conta o número de arquivos .fast5 do diretório fast5_dir
ls $1/*.fast5 | wc -l
# Conta o total de reads do diretório fast5_dir
h5ls $1/*.fast5 | wc -l

