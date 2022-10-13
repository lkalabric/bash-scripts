#!/bin/bash

# script: blast_report.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# objetivo: relatório dos resultados da calssificação taxonômica pelo blastn
# criação: 25 AGO 2021
# ultima atualização: 14 JUN 2022
# atualização: revisão do script

# Valiação da entrada de dados na linha de comando
if [ ! -f $1 ]; then
	echo "Arquivo ou diretório inválido"
	echo "Sintáxe: ./blast_report.sh <caminho_completo/barcodeXX.blastn>"
	echo "Exemplo: ./blast_report.sh ngs-library/DENV_FTA_1_fast/wf2/READSLEVEL/BLASTN/barcode01.blastn"
	exit 1
fi
FILENAME=$1

# Output filename
OUTPUTFILENAME=${FILENAME/.blastn/.blastnreport}

# Lista os vírus pesquisados
WIMPDIR="${HOME}/data/WIMP"
# echo "Lista de vírus pesquisados"
# ls ${WIMPDIR}

echo "All reads: $(wc -l ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Número de reads por taxon"  | tee -a ${OUTPUTFILENAME} 
echo "Flaviviridae"  | tee -a ${OUTPUTFILENAME} 
printf '	%s		%s\n' "DENV1	" "$(grep -c -F -f ${WIMPDIR}/denv1.wimp ${FILENAME})"  | tee -a ${OUTPUTFILENAME} 
printf '	%s		%s\n' "DENV2	" "$(grep -c -F -f ${WIMPDIR}/denv2.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "DENV3	" "$(grep -c -F -f ${WIMPDIR}/denv3.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "DENV4	" "$(grep -c -F -f ${WIMPDIR}/denv4.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "ZIKV	" "$(grep -c -F -f ${WIMPDIR}/zikv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "WNV	" "$(grep -c -F -f ${WIMPDIR}/wnv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "YFV	" "$(grep -c -F -f ${WIMPDIR}/yfv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "SLEV	" "$(grep -c -F -f ${WIMPDIR}/slev.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "JEV	" "$(grep -c -F -f ${WIMPDIR}/jev.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "ILHV	" "$(grep -c -F -f ${WIMPDIR}/ilhv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "ROCV	" "$(grep -c -F -f ${WIMPDIR}/rocv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "HCV	" "$(grep -c -F -f ${WIMPDIR}/hcv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "GB-CV	" "$(grep -c -F -f ${WIMPDIR}/gb-cv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Picornaviridae" | tee -a ${OUTPUTFILENAME} 
printf "	%s     	 	%s\n" "CV	" "$(grep -c -f ${WIMPDIR}/cv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "EV-A	" "$(grep -c -F -f ${WIMPDIR}/ev-a.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "EV-B	" "$(grep -c -F -f ${WIMPDIR}/ev-b.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "EV-C/PV	" "$(grep -c -F -f ${WIMPDIR}/ev-c.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "EV-D	" "$(grep -c -F -f ${WIMPDIR}/ev-d.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "EV-E2H	" "$(grep -c -F -f ${WIMPDIR}/ev-e2h.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "RV	" "$(grep -c -F -f ${WIMPDIR}/rv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "HRV	" "$(grep -c -F -f ${WIMPDIR}/hrv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "PV	" "$(grep -c -F -f ${WIMPDIR}/pv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s	%s\n" "OutrosEV	" "$(grep -c -F -f ${WIMPDIR}/outrosev.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Herpexviridae" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "HHV	" "$(grep -c -F -f ${WIMPDIR}/hhv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Hepeviridae" | tee -a ${OUTPUTFILENAME} 
printf "	%s	%s\n" "HEV-group	" "$(grep -c -F -f ${WIMPDIR}/hev.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Caliciviridae" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "NoV	" "$(grep -c -F -f ${WIMPDIR}/nov.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Togaviridae" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "CHIKV" "$(grep -c -F -f ${WIMPDIR}/chikv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "MAYV	" "$(grep -c -F -f ${WIMPDIR}/mayv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "SINV	" "$(grep -c -F -f ${WIMPDIR}/sinv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Paramyxoviridae" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "MUV	" "$(grep -c -F -f ${WIMPDIR}/muv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Coronaviridae" | tee -a ${OUTPUTFILENAME} 
printf "	%s	%s\n" "SARS-CoV-2	" "$(grep -c -F -f ${WIMPDIR}/sars-cov-2.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Retroviridae" | tee -a ${OUTPUTFILENAME} 
printf "	%s      		%s\n" "HIV	" "$(grep -c -F -f ${WIMPDIR}/hiv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
printf "	%s		%s\n" "HTLV	" "$(grep -c -F -f ${WIMPDIR}/htlv.wimp ${OUTPUTFILENAME})" | tee -a ${OUTPUTFILENAME} 
echo "Hepadnaviridae" | tee -a ${OUTPUTFILENAME} 
printf "	%s		%s\n" "HBV	" "$(grep -c -F -f ${WIMPDIR}/hbv.wimp ${FILENAME})" | tee -a ${OUTPUTFILENAME} 
