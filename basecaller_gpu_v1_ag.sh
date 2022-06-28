#!/bin/bash

# script: basecaller_gpu_v1_ag.sh
# author: Alessandra Gonzalez <le.sandragonzalez@gmail.com>
# institution: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# date: 18 AGO 2021
# útlima atualização: 28 JUN 2022

# Nome dos dados passado na linha de comando
RUNNAME=$1 
MODELO=$2

if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados e ou o modelo de análise!"
	exit 0
fi
# Caminho de INPUT dos dados fast5
 INPUT_DIR="$HOME/data/$RUNNAME"
 
 if [ ! -d $INPUT_DIR ]; then
 	echo "Pasta de dados não encontrada!"
 	exit 0
 fi
# Caminho de OUTPUT dos resultados da análise NGS
mkdir -p "${HOME}/ngs-analysis/${RUNNAME}_${MODELO}"
SAVE_DIR="${HOME}/ngs-analysis/${RUNNAME}_${MODELO}"
# Flow cell and kit OR config file
# FLOWCELL="FLO-MIN106"
# KIT="SQK-LSK109"
# Arquivo _fast para benckmark apenas e _hac para as análises
# CONFIG="dna_r9.4.1_450bps_fast.cfg" # "dna_r9.4.1_450bps_hac.cfg"
CONFIG="dna_r9.4.1_450bps_${MODELO}.cfg"

# Parametros para otimização da GPU
# --gpu_runners_per_device, --chunk_size, --chunks_per_runner, -x

# Comando para guppy_basecaller usando GPU do LAPTOP-Yale e parâmetros otimizados a partir do benchmark 
guppy_basecaller -r -i ${INPUT_DIR} -s ${SAVE_DIR} -c ${CONFIG} -x auto --gpu_runners_per_device 4 --chunk_size 1000 --chunks_per_runner 50 --verbose_logs
