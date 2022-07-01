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

barcode_bm () {
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
}

barcode_bm
exit

# Demultiplex Benchmarks
#bm1 
echo -e "\nExecutando barcoder_bm1..."
[ -d "${DEMUXDIR}/bm1" ] && rm -r "${DEMUXDIR}/bm1"; mkdir -vp "${DEMUXDIR}/bm1"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm1" --arrangements_files ${ARRANGEMENTS} --trim_barcodes
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm1/barcoding_summary.txt" -o "${DEMUXDIR}/${RUNNAME}_bm1_pycoqc.html" --report_title "${RUNNAME}_bm1" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm2 
echo -e "\nExecutando barcoder_bm2..."
[ -d "${DEMUXDIR}/bm2" ] && rm -r "${DEMUXDIR}/bm2"; mkdir -vp "${DEMUXDIR}/bm2"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm2" --arrangements_files ${ARRANGEMENTS} --detect_mid_strand_adapter --trim_barcodes
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm2/barcoding_summary.txt" -o "${DEMUXDIR}/${RUNNAME}_bm2_pycoqc.html" --report_title "${RUNNAME}_bm2" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm3 
echo -e "\nExecutando barcoder_bm3..."
[ -d "${DEMUXDIR}/bm3" ] && rm -r "${DEMUXDIR}/bm3"; mkdir -vp "${DEMUXDIR}/bm3"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm3" --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --trim_barcodes
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm3/barcoding_summary.txt" -o "${DEMUXDIR}/${RUNNAME}_bm3_pycoqc.html" --report_title "${RUNNAME}_bm3" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm4 
echo -e "\nExecutando barcoder_bm4..."
[ -d "${DEMUXDIR}/bm4" ] && rm -r "${DEMUXDIR}/bm4"; mkdir -vp "${DEMUXDIR}/bm4"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm4" --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --detect_mid_strand_barcodes --trim_barcodes  
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm4/barcoding_summary.txt" -o "${DEMUXDIR}/${RUNNAME}_bm4_pycoqc.html" --report_title "${RUNNAME}_bm4" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm5
echo -e "\nExecutando barcoder_bm5..."
[ -d "${DEMUXDIR}/bm5" ] && rm -r "${DEMUXDIR}/bm5"; mkdir -vp "${DEMUXDIR}/bm5"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm5" --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --detect_mid_strand_adapter --trim_barcodes  
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm5/barcoding_summary.txt" -o "${DEMUXDIR}/${RUNNAME}_bm5_pycoqc.html" --report_title "${RUNNAME}_bm5" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm6 
echo -e "\nExecutando barcoder_bm6..."
[ -d "${DEMUXDIR}/bm6" ] && rm -r "${DEMUXDIR}/bm6"; mkdir -vp "${DEMUXDIR}/bm6"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm6" --barcode_kits EXP-NBD104 --require_barcodes_both_ends --detect_mid_strand_barcodes --trim_barcodes  
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm6/barcoding_summary.txt" -o "${DEMUXDIR}/${RUNNAME}_bm6_pycoqc.html" --report_title "${RUNNAME}_bm6" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

exit
