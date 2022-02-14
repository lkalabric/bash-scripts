#!/bin/bash

# author: Alessandra Gonzalez <le.sandragonzalez@gmail.com>
# institution: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# date: 18 AGO 2021

# Nome dos dados passado pelo teclado
RUNNAME=$1
[ ! f "${RUNNAME}" ] && RUNNAME="DENV_FTA_1" 
echo "Análise da biblioteca ${RUNNAME}!"

# Caminho de INPUT dos dados fast5
INPUT_DIR="$HOME/data/$RUNNAME"
if [ ! -d $INPUT_DIR ]; then
	echo "Pasta de dados não encontrada!"
 	exit 0
fi

# Caminho de OUTPUT dos resultados da análise NGS
echo "Preparando as pastas para (re-)análise dos dados..."
SAVE_DIR="${HOME}/ngs-analysis/${RUNNAME}_bm"

# Remove os dados da última análise
rm -r ${SAVE_DIR}
[ ! -d "${SAVE_DIR}" ] && mkdir -vp ${SAVE_DIR}

# Flow cell and kit OR config file
# FLOWCELL="FLO-MIN106"
# KIT="SQK-LSK109"
# Arquivo _fast para benckmark apenas e _hac para as análises
CONFIG="dna_r9.4.1_450bps_hac.cfg" #"dna_r9.4.1_450bps_fast.cfg"

# Parametros para otimização da CPU
# --num_callers, --num_threads_per_caller

# Parametros para otimização da GPU
# --gpu_runners_per_device, --chunk_size, --chunks_per_runner, -x

# Comando para guppy_basecaller usando GPU e parâmetros padrões do arquivo de configuração 
#guppy_basecaller -r -i ${INPUT_DIR} -s ${SAVE_DIR} -c ${CONFIG} -x auto --verbose_logs

# Parâmtros estimados por Laise
# GPUPERDEVICE="$(nvidia-settings -q TotalDedicatedGPUMemory -t | awk '{print $1/10/256}' | awk -F, '{print $1}')"
#THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Comando para guppy_basecaller usando GPU e parâmetros otimizados por Laise_hac
#guppy_basecaller -r -i ${INPUT_DIR} -s ${SAVE_DIR} -c ${CONFIG} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --num_callers ${THREADS} --verbose_logs

# Parametros para benchmark (bm)
GPUPERDEVICE=8
CHUNKSIZE=2000
CHUNKPERRUNNER=256
THREADS=4

# Comando para guppy_basecaller usando GPU e configuração do benchmark
guppy_basecaller -r -i ${INPUT_DIR} -s ${SAVE_DIR} -c ${CONFIG} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --chunk_size ${CHUNKSIZE} --chunks_per_runner ${CHUNKPERRUNNER} --num_callers ${THREADS} --verbose_logs
