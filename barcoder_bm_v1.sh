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
WF=$3		# Workflow de bioinformatica 1, 2 ou 31

if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, número do worflow ou modelo Guppy Basecaller!"
	echo "Sintáxe: ./metagenomic_v31_lk.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,31>"
	exit 0
fi
# Caminho de INPUT dos dados fast5
RAWDIR="${HOME}/data/${RUNNAME}" # Se análise começar com o Barcoder
if [ ! -d $RAWDIR ]; then
	echo "Pasta de dados não encontrada!"
	exit 0
fi

# Caminho de INPUT dos bancos de dados
REFSEQDIR=${HOME}/data/REFSEQ
HUMANREFDIR=${HOME}/data/GRCh38
BLASTDBDIR="${HOME}/data/BLAST_DB"
KRAKENDB="${HOME}/data/KRAKEN2_DB" # Substituir pelo nosso banco de dados se necessário KRAKEN2_USER_DB

# Caminhos de OUTPUT das análises
echo "Preparando pastas para (re-)análise dos dados..."
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
# rm -r ${RESULTSDIR}
[ ! -d "${RESULTSDIR}" ] && mkdir -vp ${RESULTSDIR}
BASECALLDIR="${RESULTSDIR}/BASECALL"
DEMUXDIR="${RESULTSDIR}/DEMUX"
DEMUXCATDIR="${RESULTSDIR}/DEMUX_CAT"
CUTADAPTDIR="${RESULTSDIR}/wf${WF}/CUTADAPT"
NANOFILTDIR="${RESULTSDIR}/wf${WF}/NANOFILT"
PRINSEQDIR="${RESULTSDIR}/wf${WF}/PRINSEQ"
QUERYDIR="${RESULTSDIR}/wf${WF}/QUERY"
BLASTDIR="${RESULTSDIR}/wf${WF}/BLAST"
READSLEVELDIR="${RESULTSDIR}/wf${WF}/READS_LEVEL"
CONTIGLEVELDIR="${RESULTSDIR}/wf${WF}/CONTIGS_LEVEL"
ASSEMBLYDIR="${RESULTSDIR}/wf${WF}/ASSEMBLY"

# Parâmetros Guppy basecaller (ONT)
CONFIG="dna_r9.4.1_450bps_${MODEL}.cfg" #dna_r9.4.1_450bps_fast.cfg dna_r9.4.1_450bps_hac.cfg dna_r9.4.1_450bps_sup.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
 
# Parâmetros para otimização do Guppy basecaller para modelo fast utilizadando o LAPTOP-Yale (benckmark)
GPUPERDEVICE=4		
CHUNCKSIZE=1000		
CHUNKPERRUNNER=50	

# Parâmetro de otimização minimap2, samtools, racon e kraken2
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros de qualidade mínima
QSCORE=9
LENGTH=100

# Parâmetros Barcoder ou Cutadapt para remoção do primer
TRIMADAPTER=18
PRIMER="GTTTCCCACTGGAGGATA"

# Parâmetros minimap2 
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -P ${HUMANREFDIR}
HUMANREFSEQ="${HUMANREFDIR}/GRCh38.p13.genome.fa.gz"
HUMANREFMMI="${HUMANREFDIR}/GRCh38.p13.genome.mmi"
# Cria o arquivo índice do genoma humano para reduzir o tempo de alinhamento
if [ ! -f $HUMANREFMMI ]; then
	minimap2 -d $HUMANREFMMI $HUMANREFSEQ
fi

# if false; then # Desvio para execução rápida


# Step 2 - Demultiplex & adapter removal
#bm1 
echo -e "\nExecutando barcoder_bm1..."
[ -d "$DEMUXDIR/bm1" ] && rm -r "$DEMUXDIR/bm1"; mkdir -vp "$DEMUXDIR/bm1"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm1" --arrangements_files ${ARRANGEMENTS} --trim_barcodes
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm1/barcoding_summary.txt" -o "${RESULTSDIR}/${RUNNAME}_bm1_pycoqc.html" --report_title "${RUNNAME}_bm1" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm2 
echo -e "\nExecutando barcoder_bm2..."
[ -d "$DEMUXDIR/bm2" ] && rm -r "$DEMUXDIR/bm2"; mkdir -vp "$DEMUXDIR/bm2"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm2" --arrangements_files ${ARRANGEMENTS} --detect_mid_strand_adapter --trim_barcodes
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm2/barcoding_summary.txt" -o "${RESULTSDIR}/${RUNNAME}_bm2_pycoqc.html" --report_title "${RUNNAME}_bm2" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm3 
echo -e "\nExecutando barcoder_bm3..."
[ -d "$DEMUXDIR/bm3" ] && rm -r "$DEMUXDIR/bm3"; mkdir -vp "$DEMUXDIR/bm3"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm3" --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --trim_barcodes
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm3/barcoding_summary.txt" -o "${RESULTSDIR}/${RUNNAME}_bm3_pycoqc.html" --report_title "${RUNNAME}_bm3" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm4 
echo -e "\nExecutando barcoder_bm4..."
[ -d "$DEMUXDIR/bm4" ] && rm -r "$DEMUXDIR/bm4"; mkdir -vp "$DEMUXDIR/bm4"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm4" --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --detect_mid_strand_barcodes --trim_barcodes  
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm4/barcoding_summary.txt" -o "${RESULTSDIR}/${RUNNAME}_bm4_pycoqc.html" --report_title "${RUNNAME}_bm4" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm5
echo -e "\nExecutando barcoder_bm5..."
[ -d "$DEMUXDIR/bm5" ] && rm -r "$DEMUXDIR/bm5"; mkdir -vp "$DEMUXDIR/bm5"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm5" --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --detect_mid_strand_adapter --trim_barcodes  
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm5/barcoding_summary.txt" -o "${RESULTSDIR}/${RUNNAME}_bm5_pycoqc.html" --report_title "${RUNNAME}_bm5" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

#bm6 
echo -e "\nExecutando barcoder_bm6..."
[ -d "$DEMUXDIR/bm6" ] && rm -r "$DEMUXDIR/bm6"; mkdir -vp "$DEMUXDIR/bm6"
guppy_barcoder -r -i "${BASECALLDIR}/pass" -s "${DEMUXDIR}/bm6" --barcode_kits EXP-NBD104 --require_barcodes_both_ends --detect_mid_strand_barcodes --trim_barcodes  
pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/bm6/barcoding_summary.txt" -o "${RESULTSDIR}/${RUNNAME}_bm6_pycoqc.html" --report_title "${RUNNAME}_bm6" --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}

exit

