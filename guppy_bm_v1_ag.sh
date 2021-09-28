#!/bin/bash

# author: Alessandra Gonzalez <le.sandragonzalez@gmail.com>
# institution: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# date: 18 AGO 2021

# Nome dos dados passado pelo teclado
RUNNAME=$1 
if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados!"
	exit 0
fi
# Caminho de INPUT dos dados fast5
 INPUT_DIR="$HOME/data/$RUNNAME"
 
 if [ ! -d $INPUT_DIR ]; then
 	echo "Pasta de dados não encontrada!"
 	exit 0
 fi
# Caminho de OUTPUT dos resultados da análise NGS
mkdir -p "${HOME}/ngs-analysis/${RUNNAME}_bm"
SAVE_DIR="${HOME}/ngs-analysis/${RUNNAME}_bm"
# Flow cell and kit OR config file
# FLOWCELL="FLO-MIN106"
# KIT="SQK-LSK109"
# Arquivo _fast para benckmark apenas e _hac para as análises
CONFIG="dna_r9.4.1_450bps_hac.cfg" # "dna_r9.4.1_450bps_hac.cfg"

# Parametros para otimização da CPU
# --num_callers, --num_threads_per_caller

# Parametros para otimização da GPU
# --gpu_runners_per_device, --chunk_size, --chunks_per_runner, -x

# Comando para guppy_basecaller usando GPU e parâmetros padrões do arquivo de configuração 
guppy_basecaller -r -i ${INPUT_DIR} -s ${SAVE_DIR} -c ${CONFIG} -x auto --verbose_logs

# Comando para guppy_basecaller usando GPU e parâmetros otimizados a partir do benchmark 
#guppy_basecaller -r -i ${INPUT_DIR} -s ${SAVE_DIR} -c ${CONFIG} -x auto --gpu_runners_per_device 4 --chunk_size 1000 --chunks_per_runner 50 --verbose_logs

# Comando para guppy_basecaller usando GPU e parâmetros otimizados a partir do benchmark 
#guppy_basecaller -r -i ${INPUT_DIR} -s ${SAVE_DIR} -c ${CONFIG} -x auto --gpu_runners_per_device 4 --chunk_size 1000 --chunks_per_runner 50 --num_callers 12 --verbose_logs



