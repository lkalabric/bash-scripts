#!/bin/bash

# script: blast_report.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# objetivo: relatório dos resultados da calssificação taxonômica pelo blastn
# criação: 25 AGO 2021
# ultima atualização: 14 JUN 2022
# atualização: revisão do script

# Valiação da entrada de dados na linha de comando
die () {
    echo >&2 "$@"
    exit 0
}

case $# in
	1)
		if [ ! -f $1 ]; then
			echo "Arquivo ou diretório inválido"
			echo "Sintáxe: ./blast_report.sh <caminho_completo/barcode01.blastn>"
			echo "Exemplo: ./blast_report.sh ngs-library/DENV_FTA_1_fast/wf2/BLAST/barcode01.blastn"
			exit 1
		fi
		$FILENAME=$1
	;;
	2)
		RUNNAME=$1
		BARCODE=$2
		RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}"
		BLASTDIR="${RESULTSDIR}/wf2/BLAST"
		echo "Lista de taxons no BLAST_DB"
		FILENAME="${BLASTDIR}/${BARCODE}.blastn"
		if [ ! -f ${FILENAME} ]; then	
			echo "Falta o nome da biblioteca_model e/ou do barcodeXX"
			echo "Sintáxe: ./blast_report.sh <BIBLIOTECA_MODEL> <BARCODE>"
			echo "Exemplo: ./blast_report.sh DENV_FTA_1_fast barcode01"
			exit 2
		fi
	;;
	*)
		die "Mínimo de 1 e máximo de 2 argumentos são requiridos, $# provido"
	;;
esac

# Lista os vírus pesquisados
WIMPDIR="${HOME}/data/WIMP"
echo "Lista de vírus pesquisados"
ls ${WIMPDIR}

echo "All reads: $(wc -l ${FILENAME})" 
echo "Número de reads por taxon"
echo "Flaviviridae"
printf '	%s		%s\n' "DENV1	" "$(grep -c -F -f ${WIMPDIR}/denv1.wimp ${FILENAME})"
printf '	%s		%s\n' "DENV2	" "$(grep -c -F -f ${WIMPDIR}/denv2.wimp ${FILENAME})"
printf "	%s      		%s\n" "DENV3	" "$(grep -c -F -f ${WIMPDIR}/denv3.wimp ${FILENAME})"
printf "	%s      		%s\n" "DENV4	" "$(grep -c -F -f ${WIMPDIR}/denv4.wimp ${FILENAME})"
printf "	%s      		%s\n" "ZIKV	" "$(grep -c -F -f ${WIMPDIR}/zikv.wimp ${FILENAME})"
printf "	%s      		%s\n" "WNV	" "$(grep -c -F -f ${WIMPDIR}/wnv.wimp ${FILENAME})"
printf "	%s      		%s\n" "YFV	" "$(grep -c -F -f ${WIMPDIR}/yfv.wimp ${FILENAME})"
printf "	%s      		%s\n" "SLEV	" "$(grep -c -F -f ${WIMPDIR}/slev.wimp ${FILENAME})"
printf "	%s      		%s\n" "JEV	" "$(grep -c -F -f ${WIMPDIR}/jev.wimp ${FILENAME})"
printf "	%s      		%s\n" "ILHV	" "$(grep -c -F -f ${WIMPDIR}/ilhv.wimp ${FILENAME})"
printf "	%s      		%s\n" "ROCV	" "$(grep -c -F -f ${WIMPDIR}/rocv.wimp ${FILENAME})"
printf "	%s      		%s\n" "HCV	" "$(grep -c -F -f ${WIMPDIR}/hcv.wimp ${FILENAME})"
echo "Picornaviridae"
printf "	%s     	 	%s\n" "CV	" "$(grep -c -f ${WIMPDIR}/cv.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-A	" "$(grep -c -F -f ${WIMPDIR}/ev-a.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-B	" "$(grep -c -F -f ${WIMPDIR}/ev-b.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-C/PV	" "$(grep -c -F -f ${WIMPDIR}/ev-c.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-D	" "$(grep -c -F -f ${WIMPDIR}/ev-d.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-E2H	" "$(grep -c -F -f ${WIMPDIR}/ev-e2h.wimp ${FILENAME})"
printf "	%s      		%s\n" "RV	" "$(grep -c -F -f ${WIMPDIR}/rv.wimp ${FILENAME})"
printf "	%s      		%s\n" "HRV	" "$(grep -c -F -f ${WIMPDIR}/hrv.wimp ${FILENAME})"
printf "	%s      		%s\n" "PV	" "$(grep -c -F -f ${WIMPDIR}/pv.wimp ${FILENAME})"
printf "	%s	%s\n" "OutrosEV	" "$(grep -c -F -f ${WIMPDIR}/outrosev.wimp ${FILENAME})"
echo "Herpexviridae"
printf "	%s      		%s\n" "HHV	" "$(grep -c -F -f ${WIMPDIR}/hhv.wimp ${FILENAME})"
echo "Hepeviridae"
printf "	%s	%s\n" "HEV-group	" "$(grep -c -F -f ${WIMPDIR}/hev.wimp ${FILENAME})"
echo "Caliciviridae"
printf "	%s      		%s\n" "NoV	" "$(grep -c -F -f ${WIMPDIR}/nov.wimp ${FILENAME})"
echo "Togaviridae"
printf "	%s      		%s\n" "CHIKV" "$(grep -c -F -f ${WIMPDIR}/chikv.wimp ${FILENAME})"
printf "	%s      		%s\n" "MAYV	" "$(grep -c -F -f ${WIMPDIR}/mayv.wimp ${FILENAME})"
printf "	%s      		%s\n" "SINV	" "$(grep -c -F -f ${WIMPDIR}/sinv.wimp ${FILENAME})"
echo "Paramyxoviridae"
printf "	%s      		%s\n" "MUV	" "$(grep -c -F -f ${WIMPDIR}/muv.wimp ${FILENAME})"
echo "Coronaviridae"
printf "	%s	%s\n" "SARS-CoV-2	" "$(grep -c -F -f ${WIMPDIR}/sars-cov-2.wimp ${FILENAME})"
echo "Retroviridae"
printf "	%s      		%s\n" "HIV	" "$(grep -c -F -f ${WIMPDIR}/hiv.wimp ${FILENAME})"
printf "	%s		%s\n" "HTLV	" "$(grep -c -F -f ${WIMPDIR}/htlv.wimp ${FILENAME})"
echo "Hepadnaviridae"
printf "	%s		%s\n" "HBV	" "$(grep -c -F -f ${WIMPDIR}/hbv.wimp ${FILENAME})"
