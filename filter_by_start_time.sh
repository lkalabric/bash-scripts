#!/bin/bash

# script: filter_by_star_time.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 13 ABR 2022
# última atualização: 13 ABR 2022
# versão 0: 
# link: https://www.biostars.org/p/389679/

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
STAR_DATE=#3	# Data da corrida

if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, número do worflow ou modelo Guppy Basecaller!"
	echo "Sintáxe: ./filter_by_start_time.sh <LIBRARY> <MODELO:fast,hac,sup>"
	exit 0
fi

# Caminhos de OUTPUT das análises
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
FILTER_BY_START_TIMEDIR="${RESULTSDIR}/FILTER_BY_START_TIME"
[ ! -d "${FILTER_BY_START_TIMEDIR}" ] && mkdir -vp "${FILTER_BY_START_TIMEDIR}"

# Filter_by_start_time
# Tempos para análise
declare -a START_TIME=(01 02 04 12 24 72)
for i in $(find ${DEMUXCATDIR} -type f -exec basename {} .fastq \;); do
	echo -e "\nContando as reads do arquivo ${DEMUXCATDIR}/${i}.fastq..."
	wc -l "${DEMUXCATDIR}/${i}.fastq"
	for j in "${START_TIME[@]}"; do
		echo -e "\nExecutando filter_by_start_time ${j}..."
		echo -e "${START_DATE}T${j}"
		echo -e "${DEMUXCATDIR}/${i}.fastq"
		echo -e "${FILTER_BY_START_TIMEDIR}/${i}.${j}.fastq"
		read -p "Pressione qualquer tecla para continuar..."
		grep -A3 "${START_DATE}T${j}" "${DEMUXCATDIR}/${i}.fastq" > "${FILTER_BY_START_TIMEDIR}/${i}.${j}.fastq"
		wc -l "${FILTER_BY_START_TIMEDIR}/${i}.${j}.fastq"
	done
done

