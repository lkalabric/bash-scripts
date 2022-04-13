#!/bin/bash

# script: filter_by_star_time.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 13 ABR 2022
# última atualização: 13 ABR 2022
# versão 0: 

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=$3		# Workflow de bioinformatica 1, 2 ou 31

if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, número do worflow ou modelo Guppy Basecaller!"
	echo "Sintáxe: ./metagenomic_v31_lk.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,31>"
	exit 0
fi

# Caminhos de OUTPUT das análises
echo "Preparando pastas para (re-)análise dos dados..."
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
[ ! -d "${RESULTSDIR}/FILTER_BY_START_TIME" ] && mkdir -vp "${RESULTSDIR}/FILTER_BY_START_TIME"
BASECALLDIR="${RESULTSDIR}/BASECALL"
DEMUXDIR="${RESULTSDIR}/DEMUX"
DEMUXCATDIR="${RESULTSDIR}/DEMUX_CAT"
FILTER_BY_START_TIMEDIR="${RESULTSDIR}/FILTER_BY_START_TIME"
[ ! -d "${FILTER_BY_START_TIMEDIR}" ] && mkdir -vp "${FILTER_BY_START_TIMEDIR}"

# Parâmetro de otimização minimap2, samtools, racon e kraken2
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros de qualidade mínima
QSCORE=9
LENGTH=100

# Filter_by_start_time
# Tempos para análise
declare -a START_TIME=(1 2 4 12 24 72)
echo -e "\nExecutando cutadapt..."
[ ! -d ${CUTADAPTDIR} ] && mkdir -vp ${CUTADAPTDIR}
for i in $(find ${DEMUXCATDIR} -type f -exec basename {} .fastq \;); do
  wc -l "${DEMUXCATDIR}/${i}.fastq"
  for j in $({START_TIME[@]}) do
    grep -A3 "2021-06-02T${START_TIME}" "${DEMUXCATDIR}/${i}.fastq" > "${FILTER_BY_START_TIMEDIR}/${i}.${j}.fastq"
    wc -l "${FILTER_BY_START_TIMEDIR}/${i}.${j}.fastq"
  done
 done



grep -A3 "2019-03-11T16" x.fastq > test.fastq
