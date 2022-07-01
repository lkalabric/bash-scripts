#!/bin/bash

# script: barcoder_bm_v1.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 25 AGO 2021
# última atualização: 30 JUN 2022
# versão 1: separou a execução do barcoder do script metagenomic_v31_lk.sh e permitiu executar as diferentes versões de barcoder de uma só vez!

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup

if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, número do worflow ou modelo Guppy Basecaller!"
	echo "Sintáxe: ./barcoder_bm_v1.sh <LIBRARY> <MODELO:fast,hac,sup>"
	exit 0
fi

# Caminhos de OUTPUT das análises
echo "Preparando pastas para (re-)análise dos dados..."
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
# rm -r ${RESULTSDIR}
[ ! -d "${RESULTSDIR}" ] && mkdir -vp ${RESULTSDIR}
BASECALLDIR="${RESULTSDIR}/BASECALL"
DEMUXDIR="${RESULTSDIR}/DEMUX_bm"

# Parâmetros Guppy basecaller (ONT)
CONFIG="dna_r9.4.1_450bps_${MODEL}.cfg" #dna_r9.4.1_450bps_fast.cfg dna_r9.4.1_450bps_hac.cfg dna_r9.4.1_450bps_sup.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
 
# Parâmetros de qualidade mínima
QSCORE=9
LENGTH=100

# if false; then # Desvio para execução rápida

BARCODE_ARGS=(
	'--trim_barcodes'
	'--detect_mid_strand_adapter --trim_barcodes'
	'--require_barcodes_both_ends --trim_barcodes'
	'--require_barcodes_both_ends --detect_mid_strand_barcodes --trim_barcodes'
	'--require_barcodes_both_ends  --detect_mid_strand_barcodes --detect_mid_strand_adapter --trim_barcodes  '
	'--barcode_kits EXP-NBD104 --require_barcodes_both_ends --detect_mid_strand_barcodes --trim_barcodes'
)
# Contador do número de BMs
BM=1
for i in "${BARCODE_ARGS[@]}"; do
	echo -e "\nExecutando barcoder_bm${BM}..."
	[ -d "${DEMUXDIR}/bm${BM}" ] && rm -r "${DEMUXDIR}/bm${BM}"; mkdir -vp "${DEMUXDIR}/bm${BM}"
	guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm${BM}" --arrangements_files ${ARRANGEMENTS} ${i}
	pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm${BM}/barcoding_summary.txt" -o "${DEMUXDIR}/${RUNNAME}_bm${BM}_pycoqc.html" --report_title "${RUNNAME}_bm1" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}
	BM=$((BM+1))
done
