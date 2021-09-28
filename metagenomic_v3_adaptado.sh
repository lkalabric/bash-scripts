#!/bin/bash

# author: Laise de Moraes <laisepaixao@live.com>
# adaptado: Luciano Kalabric & Alessandra Gonzalez
# institution: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# URL: https://lpmor22.github.io
# date: 12 JUN 2021

# Nome dos dados passado pelo teclado
RUNNAME=$1 
if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados!"
	exit 0
fi
# Caminho de INPUT dos dados fast5
RAWDIR="$HOME/ngs-analysis/$RUNNAME"
 
if [ ! -d $RAWDIR ]; then
	echo "Pasta de dados não encontrada!"
 	exit 0
fi
CONFIG="dna_r9.4.1_450bps_hac.cfg" #dna_r9.4.1_450bps_fast.cfg dna_r9.4.1_450bps_sup.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
TRIMADAPTER="18"

HUMANREFSEQ="$HOME/data/GRCh38/GRCh38.p13.genome.fa.gz"
HUMANREFMMI="$HOME/data/GRCh38/GRCh38.p13.genome.mmi"
KRAKENDB="$HOME/data/MINIKRAKEN2_DB" # Substituir pelo nosso banco de dados se necessário KRAKEN2_USER_DB

GPUPERDEVICE="$(nvidia-settings -q TotalDedicatedGPUMemory -t | awk '{print $1/10/256}' | awk -F, '{print $1}')"
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"
BASECALLDIR="${RAWDIR}/BASECALL"
DEMUXDIR="${RAWDIR}/DEMUX"
DEMUXCATDIR="${RAWDIR}/DEMUX_CAT"
READLEVELDIR="${RAWDIR}/LEVEL_READS"
CONTIGLEVELDIR="${RAWDIR}/LEVEL_CONTIGS"
ASSEMBLYDIR="${RAWDIR}/ASSEMBLY"


if false; then # Desvio para debug1
# Esta etapa está sendo realizada à parte pelo script guppy_bm_v1_ag.sh

# Step 1 - Basecalling
echo Executando guppy_basecaller...
# guppy_basecaller -r -i ${RAWDIR} -s ${BASECALLDIR} -c ${CONFIG} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --num_callers ${THREADS} --verbose_logs
guppy_basecaller -r -i ${RAWDIR} -s ${BASECALLDIR} -c ${CONFIG} -x auto  --gpu_runners_per_device 4 --chunk_size 1000 --chunks_per_runner 50

# Step 2 - Demultiplex & adapter removal
echo "Executando guppy_barcoder..."
mkdir -p "${RAWDIR}/DEMUX"
# guppy_barcoder -r -i ${BASECALLDIR} -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER} -t ${THREADS} -x auto
guppy_barcoder -r -i ${BASECALLDIR} -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER} -x auto

fi # Desvio para debub 1
#exit 0 # Fim do debub 1

# Step 3 - QC
echo "Executando pycoQC..."
source activate pycoqc
# Com filtro QC
pycoQC --min_pass_qual 9 --min_pass_len 100 -q -f ${BASECALLDIR}/sequencing_summary.txt -b ${DEMUXDIR}/barcoding_summary.txt -o ${RAWDIR}/${RUNNAME}_qc.html --report_title ${RUNNAME}

# Concatena todos arquivos .fastq de cada barcode em um arquivo .fastq único
[ ! -d "${DEMUXCATDIR}" ] && mkdir -vp ${DEMUXCATDIR}
for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
    [ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > ${DEMUXCATDIR}/${i}.fastq
done

# Step 4 - Remoção das reads do genoma humano
echo "Executando minimap2, samtools & racon..."
[ ! -d "${READLEVELDIR}" ] && mkdir -vp ${READLEVELDIR}
[ ! -d "${ASSEMBLYDIR}" ] && mkdir -vp ${ASSEMBLYDIR}
# Cria o arquivo índice do genoma humano para reduzir o tempo de alinhamento
minimap2 -d ${HUMANREFMMI} ${HUMANREFSEQ}

# Loop para analisar todos barcodes, um de cada vez
for i in $(find ${DEMUXCATDIR} -type f -name "*.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	source activate minimap2
    	# Alinha as reads contra o arquivo indice do genoma humano e ordena os segmentos
    	minimap2 -ax map-ont -t ${THREADS} ${HUMANREFMMI} ${DEMUXCATDIR}/${i}.fastq | samtools sort -@ 12 -o ${READLEVELDIR}/${i}.sorted.bam -
    	# Indexa o arquivo para acesso mais rápido
    	samtools index -@ ${THREADS} ${READLEVELDIR}/${i}.sorted.bam
    	# Filtra os segmento não mapeados Flag 4 (-f 4) para um novo arquivo .bam 
    	samtools view -bS -f 4 ${READLEVELDIR}/${i}.sorted.bam > ${READLEVELDIR}/${i}.unmapped.sam -@ 12
	# Salva os dados no formato .fastq
	samtools fastq -f 4 ${READLEVELDIR}/${i}.unmapped.sam > ${READLEVELDIR}/${i}.unmapped.fastq -@ 12
	# Alinhar todas as reads com elas mesmas para produzir sequencias consenso a partir do overlap de reads
	minimap2 -ax ava-ont -t ${THREADS} ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.overlap.sam
	# Correção de erros a partir das sequencias consenso
 	source activate racon
 	racon -t ${THREADS} -f -u ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.overlap.sam ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.corrected.fasta
done

# Step 5 - Classificação taxonômica
echo "Executando o Kraken2..."
for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	# kraken2 --db ${KRAKENDB} --threads ${THREADS} --report ${READLEVELDIR}/${i}_report.txt --report-minimizer-data --output ${READLEVELDIR}/${i}_output.txt ${READLEVELDIR}/${i}.corrected.fasta
	echo "Carregando os dados ${READLEVELDIR}/${i}..."
	kraken2 --db ${KRAKENDB} --quick --threads ${THREADS} --report ${READLEVELDIR}/${i}_report.txt --output ${READLEVELDIR}/${i}_output.txt ${READLEVELDIR}/${i}.corrected.fasta

done

echo "Análise concluida com sucesso!"
exit 0

CHIKVREFSEQ=""
DENV1REFSEQ=""
DENV2REFSEQ=""
DENV3REFSEQ=""
DENV4REFSEQ=""
ZIKVREFSEQ=""

#source activate minimap2

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${CHIKVREFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.chikv.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.chikv.sorted.bam > ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${CHIKVREFSEQ} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.chikv -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV1REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv1.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv1.sorted.bam > ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV1REFSEQ} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv1 -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV2REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv2.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv2.sorted.bam > ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV2REFSEQ} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv2 -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV3REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv3.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv3.sorted.bam > ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV3REFSEQ} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv3 -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV4REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv4.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv4.sorted.bam > ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV4REFSEQ} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv4 -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${ZIKVREFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.zikv.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.zikv.sorted.bam > ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${ZIKVREFSEQ} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.zikv -n N -i ${i}
#done

#source activate plot

#for i in $(find ${ASSEMBLYDIR} -type f -name "*.sorted.mapped.bam" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	fastcov ${ASSEMBLYDIR}/barcode*.chikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_chikv.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.chikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_chikv_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.denv1.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv1.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.denv1.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv1_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.denv2.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv2.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.denv2.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv2_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.denv3.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv3.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.denv3.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv3_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.denv4.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv4.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.denv4.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv4_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.zikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_zikv.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.zikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_zikv_log.pdf
#done
