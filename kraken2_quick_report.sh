#!/bin/bash

# script: kraken2_quick_report.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# objetivo: relatório resumido dos resultados da classificação taxonômica pelo kraken2
# criação: 25 AGO 2021
# ultima atualização: 17 OUT 2021
# atualização: revisão do script
# requisito: arquivo taxin que contém a lista de taxons de interesse

# Valiação da entrada de dados na linha de comando
case $# in
	1)
		if [ ! -f $1 ]; then
			echo "Arquivo ou diretório inválido"
			echo "Sintáxe: ./kraken2_quick_report.sh <caminho_completo/barcodeXX._report.txt>"
			echo "Exemplo: ./kraken2_quick_report.sh ngs-library/DENV_FTA_1_fast/wf3/READS_LEVEL/barcode01_report.txt"
			exit 1
		fi
		FILENAME=$1
	;;
	2)
		RUNNAME=$1
		BARCODE=$2
		RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}"
		READSLEVELDIR="${RESULTSDIR}/wf3/READS_LEVEL"
		FILENAME="${READSLEVELDIR}/${BARCODE}_report.txt"
		if [ ! -f ${FILENAME} ]; then	
			echo "Falta o nome da biblioteca_model e/ou do barcodeXX"
			echo "Sintáxe: ./kraken2_quick_report.sh <BIBLIOTECA_MODEL> <BARCODE>"
			echo "Exemplo: ./kraken2_quick_report.sh DENV_FTA_1_hac barcode01"
			exit 2
		fi
	;;	
	*)
		echo "Mínimo de 1 e máximo de 2 argumentos são requiridos, $# provido"
		echo "Sintáxe: ./kraken2_quick_report.sh <caminho_completo/barcodeXX_report.txt>"
		echo "Sintáxe: ./kraken2_quick_report.sh <BIBLIOTECA_MODEL> <BARCODE>"
		exit 0
	;;
esac
OUTPUTFILENAME=$(echo ${FILENAME} | replace "_report.txt" "_report.kraken")

while read -r line ; do
	count=$(agrep -q -w "$line" ${FILENAME} | cut -f 2)
	echo "$line - $count" | tee -a ${OUTPUTFILENAME}
done < <(cat /home/brazil1/data/REFSEQ/taxin | tr '\t' ';')
