#!/bin/bash

# script: filter_by_star_time.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 13 ABR 2022
# última atualização: 02 DEZ 2022
# versão 2: 
# link: https://www.biostars.org/p/389679/

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=$3		# Workflow de bioinformatica 1, 2 ou 3

if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, modelo ou data da corrida!"
	echo "Sintáxe: ./filter_by_start_time.sh <LIBRARY> <MODELO:fast,hac,sup> <WF> <DATA:YYYY-MM-DD>"
	exit 0
fi

# Caminhos de INPUT das análises
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
DEMUXCATDIR="${RESULTSDIR}/wf${WF}/DEMUX_CAT"

# Caminho de OUTPUT das análises
FILTER_BY_START_TIMEDIR="${RESULTSDIR}/wf${WF}_filter/DEMUX_CAT2"
[ ! -d "${FILTER_BY_START_TIMEDIR}" ] && mkdir -vp "${FILTER_BY_START_TIMEDIR}"

# Filter_by_start_time
# Expressões regulares para filtar pelo tempo da corrida
declare -a START_TIME=(01 02 04 08 12 16 24 all)
declare -a REGEXP=("(0+[0-1])" "(0+[0-2])" "(0+[0-4])" "(0+[0-8])" "(0+[0-9]|1+[0-2])" "(0+[0-9]|1+[0-6])" "(0+[0-9]|1+[0-9]|2+[0-4])" "..")
for i in $(find ${DEMUXCATDIR} -type f -exec basename {} .fastq \; | sort); do
	echo -e "\nContando as reads do arquivo ${DEMUXCATDIR}/${i}.fastq..."
	grep "${RUNNAME}" "${DEMUXCATDIR}/${i}.fastq" | wc -l
	for ((j=0; j<=7; j++)); do
		echo -e "\nExecutando filter_by_start_time ${START_TIME[${j}]}..."
		egrep -A3 "start_time=..........T${REGEXP[${j}]}" "${DEMUXCATDIR}/${i}.fastq" > "${FILTER_BY_START_TIMEDIR}/${i}_${START_TIME[${j}]}.fastq"
		echo -e "\nContando as reads do arquivo "${FILTER_BY_START_TIMEDIR}/${i}_${START_TIME[${j}]}.fastq"..."
		grep "${RUNNAME}" "${FILTER_BY_START_TIMEDIR}/${i}_${START_TIME[${j}]}.fastq" | wc -l
	done
done

