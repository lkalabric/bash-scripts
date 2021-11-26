#!/bin/bash

# script: blast_report.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# objetivo: relatório dos resultados da calssificação taxonômica pelo blastn
# criação: 25 AGO 2021
# ultima atualização: 17 OUT 2021
# atualização: revisão do script

# Valiação da entrada de dados na linha de comando
if [[ $# -ne 2 ]]; then
	echo "Falta o nome da biblioteca_model e/ou do barcodeXX!"
	echo "Sintáxe: ./blast_report.sh <BIBLIOTECA_MODEL> <BARCODE>"
	echo "Exemplo: ./blast_report.sh DENV_FTA_1_hac barcode01"
	exit 0
fi

RUNNAME=$1
BARCODE=$2
WF=2
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}"
BLASTDIR="${RESULTSDIR}/wf${WF}/BLAST"
echo "Lista de taxons no BLAST_DB"
WIMPDIR="${HOME}/data/REFSEQ/WIMP"
ls ${WIMPDIR}
FILENAME="${BLASTDIR}/${BARCODE}.blastn"

echo "All reads: $(wc -l ${FILENAME})" 
echo "Número de reads por taxon"
echo "Flaviviridae"
printf '	%s		%s\n' "DENV1	" "$(grep -c -F -f ${WIMPDIR}/denv1.wimp ${FILENAME})"
printf '	%s		%s\n' "DENV2	" "$(grep -c -F -f ${WIMPDIR}/denv2.wimp ${FILENAME})"
printf "	%s      		%s\n" "DENV3	" "$(grep -c -F -f ${WIMPDIR}/denv3.wimp ${FILENAME})"
printf "	%s      		%s\n" "DENV4	" "$(grep -c -F -f ${WIMPDIR}/denv4.wimp ${FILENAME})"
printf "	%s      		%s\n" "ZIKV	" "$(grep -c -f ${WIMPDIR}/zikv.wimp ${FILENAME})"
printf "	%s      		%s\n" "WNV	" "$(grep -c -f ${WIMPDIR}/wnv.wimp ${FILENAME})"
printf "	%s      		%s\n" "YFV	" "$(grep -c -f ${WIMPDIR}/yfv.wimp ${FILENAME})"
printf "	%s      		%s\n" "SLEV	" "$(grep -c -f ${WIMPDIR}/slev.wimp ${FILENAME})"
printf "	%s      		%s\n" "JEV	" "$(grep -c -f ${WIMPDIR}/jev.wimp ${FILENAME})"
printf "	%s      		%s\n" "ILHV	" "$(grep -c -f ${WIMPDIR}/ilhv.wimp ${FILENAME})"
printf "	%s      		%s\n" "ROCV	" "$(grep -c -f ${WIMPDIR}/rocv.wimp ${FILENAME})"
printf "	%s      		%s\n" "HCV	" "$(grep -c -f ${WIMPDIR}/hcv.wimp ${FILENAME})"
echo "Picornaviridae"
printf "	%s     	 	%s\n" "CV	" "$(grep -c -f ${WIMPDIR}/cv.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-A	" "$(grep -c -f ${WIMPDIR}/ev-a.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-B	" "$(grep -c -f ${WIMPDIR}/ev-b.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-C/PV	" "$(grep -c -f ${WIMPDIR}/ev-c.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-D	" "$(grep -c -f ${WIMPDIR}/ev-d.wimp ${FILENAME})"
printf "	%s      		%s\n" "EV-E2H	" "$(grep -c -f ${WIMPDIR}/ev-e2h.wimp ${FILENAME})"
printf "	%s      		%s\n" "RV	" "$(grep -c -f ${WIMPDIR}/rv.wimp ${FILENAME})"
printf "	%s      		%s\n" "HRV	" "$(grep -c -f ${WIMPDIR}/hrv.wimp ${FILENAME})"
printf "	%s      		%s\n" "PV	" "$(grep -c -f ${WIMPDIR}/pv.wimp ${FILENAME})"
printf "	%s	%s\n" "OutrosEV	" "$(grep -c -f ${WIMPDIR}/outrosev.wimp ${FILENAME})"
echo "Herpexviridae"
printf "	%s      		%s\n" "HHV	" "$(grep -c -f ${WIMPDIR}/hhv.wimp ${FILENAME})"
echo "Hepeviridae"
printf "	%s	%s\n" "HEV-group	" "$(grep -c -f ${WIMPDIR}/hev.wimp ${FILENAME})"
echo "Caliciviridae"
printf "	%s      		%s\n" "NoV	" "$(grep -c -f ${WIMPDIR}/nov.wimp ${FILENAME})"
echo "Togaviridae"
printf "	%s      		%s\n" "CHIKV" "$(grep -c -f ${WIMPDIR}/chikv.wimp ${FILENAME})"
printf "	%s      		%s\n" "MAYV	" "$(grep -c -f ${WIMPDIR}/mayv.wimp ${FILENAME})"
printf "	%s      		%s\n" "SINV	" "$(grep -c -f ${WIMPDIR}/sinv.wimp ${FILENAME})"
echo "Paramyxoviridae"
printf "	%s      		%s\n" "MUV	" "$(grep -c -f ${WIMPDIR}/muv.wimp ${FILENAME})"
echo "Coronaviridae"
printf "	%s	%s\n" "SARS-CoV-2	" "$(grep -c -f ${WIMPDIR}/sars-cov-2.wimp ${FILENAME})"
echo "Retroviridade"
printf "	%s      		%s\n" "HIV	" "$(grep -c -f ${WIMPDIR}/hiv.wimp ${FILENAME})"
printf "	%s		%s\n" "HTLV	" "$(grep -c -f ${WIMPDIR}/htlv.wimp ${FILENAME})"


