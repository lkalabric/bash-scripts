#!/bin/bash

# script: coverage.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# última atualização: MAR, 17 2022
# versão 2: roda fastcov.py & samtools coverage

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=$3		#Worflow de bioinformatica

if [[ $# -eq 0 ]]; then
	echo "Falta o nome da library,modelo Guppy Basecaller ou número do workflow!"
	echo "Sintáxe: ./fastcov.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,3>"
	exit 0
fi

# Caminho INPUT dos bancos de dados
REFGENDIR="${HOME}/data/REFGEN"

# Caminhos de OUTPUT das análises
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
BASECALLDIR="${RESULTSDIR}/BASECALL"
DEMUXDIR="${RESULTSDIR}/DEMUX"
DEMUXCATDIR="${RESULTSDIR}/DEMUX_CAT"
CUTADAPTDIR="${RESULTSDIR}/wf${WF}/CUTADAPT"
NANOFILTDIR="${RESULTSDIR}/wf${WF}/NANOFILT"
PRINSEQDIR="${RESULTSDIR}/wf${WF}/PRINSEQ"
QUERYDIR="${RESULTSDIR}/wf${WF}/QUERY"
BLASTDIR="${RESULTSDIR}/wf${WF}/BLAST"
READSLEVELDIR="${RESULTSDIR}/wf${WF}/READS_LEVEL"
CONTIGSLEVELDIR="${RESULTSDIR}/wf${WF}/CONTIGS_LEVEL"
ASSEMBLYDIR="${RESULTSDIR}/wf${WF}/ASSEMBLY"
[ ! -d $ASSEMBLYDIR ] && mkdir -vp $ASSEMBLYDIR
COVERAGEDIR="${RESULTSDIR}/wf${WF}/COVERAGE"
[ ! -d $COVERAGEDIR ] && mkdir -vp $COVERAGEDIR


# Parâmetro de otimização das análises por CPU
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros para controle da qualidade mínima
QSCORE=9
LENGTH=100

# Cria o arquivo índice das sequencias referencias para mapeamento das reads pelo minimap2
for j in $(find ${REFGENDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
  echo "Gerando ídice para o genoma referencia ${j}..."
  REFGENFASTA="${REFGENDIR}/${j}.fasta"
  REFGENMMI="${REFGENDIR}/${j}.mmi"
  [ ! -f $REFGENMMI ] && minimap2 -d $REFGENMMI $REFGENFASTA
done  

# Mapeamento das sequencias em genomas referência e análise de cobertura
source activate ngs
for j in $(find ${REFGENDIR} -type f -name "*.mmi" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		echo "Mapeando ${READSLEVELDIR}/${i}.corrected.fasta..."
		minimap2 -t ${THREADS} -ax map-ont ${REFGENDIR}/${j}.mmi ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.{j}.sorted.bam -	
		samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.${j}.sorted.bam > ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam
		samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam
		samtools mpileup -A -d 0 -Q 0 ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.${j}
		samtools coverage ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam > ${COVERAGEDIR}/${i}.${j}.coverage.txt
		samtools coverage -A -w 32 ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam > ${COVERAGEDIR}/${i}.${j}.histogram.txt
		samtools depth ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam > ${COVERAGEDIR}/${i}.${j}.depth.txt
		fastcov.py ${ASSEMBLYDIR}/${i}.{j}.sorted.mapped.bam -o ${COVERAGEDIR}/${i}.{j}.fastcov.pdf -c ${COVERAGEDIR}/${i}.{j}.fastcov.txt
	done
done

exit

# Mapeamento CHIKV
#CHIKVREFGEN="${REFGENDIR}/NC_004162_CHIKV-S27.fasta"
#CHIKVREFMMI="${REFGENDIR}/NC_004162_CHIKV-S27.mmi"
#minimap2 -d $CHIKVREFMMI $CHIKVREFSEQ

#for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	echo "Mapeando ${READSLEVELDIR}/${i}.corrected.fasta..."
#	minimap2 -t ${THREADS} -ax map-ont ${CHIKVREFMMI} ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.chikv.sorted.bam -	
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.chikv.sorted.bam > ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools mpileup -A -d 0 -Q 0 ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.chikv
#done

#DENV1REFSEQ="${REFSEQDIR}/Flaviviridae/NC_001477.1_DENV1.fasta"
#DENV2REFSEQ="${REFSEQDIR}/Flaviviridae/NC_001474.2_DENV2.fasta"
#DENV3REFSEQ="${REFSEQDIR}/Flaviviridae/NC_001475.2_DENV3.fasta"
#DENV4REFSEQ="${REFSEQDIR}/Flaviviridae/NC_002640.1_DENV4.fasta"
#ZIKVREFSEQ="${REFSEQDIR}/Flaviviridae/NC_012532.1_ZIKV.fasta"

# Mapeamento DENV1
#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV1REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv1.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv1.sorted.bam > ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV1REFSEQ} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv1 -n N -i ${i}
#done

# Mapeamento DENV2
#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV2REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv2.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv2.sorted.bam > ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV2REFSEQ} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv2 -n N -i ${i}
#done

# Mapeamento DENV3
#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV3REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv3.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv3.sorted.bam > ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV3REFSEQ} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv3 -n N -i ${i}
#done

# Mapeamento DENV4
#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV4REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv4.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv4.sorted.bam > ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV4REFSEQ} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv4 -n N -i ${i}
#done

# Mapeamento ZIKV
#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${ZIKVREFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.zikv.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.zikv.sorted.bam > ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${ZIKVREFSEQ} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.zikv -n N -i ${i}
#done

source activate ngs
for i in $(find ${ASSEMBLYDIR} -type f -name "*.sorted.mapped.bam" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	echo "Determinando a cobertura do assembly ${ASSEMBLYDIR}/${i}.hcv.sorted.mapped.bam..."
	fastcov.py ${ASSEMBLYDIR}/${i}.hcv.sorted.mapped.bam -o ${COVERAGEDIR}/${i}.hcv.fastcov.pdf -c ${COVERAGEDIR}/${i}.hcv.fastcov.txt
	fastcov.py -l ${ASSEMBLYDIR}/${i}.hcv.sorted.mapped.bam -o ${COVERAGEDIR}/${i}.hcv_log.fastcov.pdf
#	fastcov.py ${ASSEMBLYDIR}/barcode*.denv1.sorted.mapped.bam -o ${COVERAGEDIR}/denv1.fastcov.pdf
#	fastcov.py -l ${ASSEMBLYDIR}/barcode*.denv1.sorted.mapped.bam -o ${COVERAGEDIR}/denv1_log.fastcov.pdf
#	fastcov.py ${ASSEMBLYDIR}/barcode*.denv2.sorted.mapped.bam -o ${COVERAGEDIR}/denv2.fastcov.pdf
#	fastcov.py -l ${ASSEMBLYDIR}/barcode*.denv2.sorted.mapped.bam -o ${COVERAGEDIR}/denv2_log.fastcov.pdf
#	fastcov.py ${ASSEMBLYDIR}/barcode*.denv3.sorted.mapped.bam -o ${COVERAGEDIR}/denv3.fastcov.pdf
#	fastcov.py -l ${ASSEMBLYDIR}/barcode*.denv3.sorted.mapped.bam -o ${COVERAGEDIR}/denv3_log.fastcov.pdf
#	fastcov.py ${ASSEMBLYDIR}/barcode*.denv4.sorted.mapped.bam -o ${COVERAGEDIR}/denv4.fastcov.pdf
#	fastcov.py -l ${ASSEMBLYDIR}/barcode*.denv4.sorted.mapped.bam -o ${COVERAGEDIR}/denv4_log.fastcov.pdf
#	fastcov.py ${ASSEMBLYDIR}/barcode*.zikv.sorted.mapped.bam -o ${COVERAGEDIR}/zikv.fastcov.pdf
#	fastcov.py -l ${ASSEMBLYDIR}/barcode*.zikv.sorted.mapped.bam -o ${COVERAGEDIR}/zikv_log.fastcov.pdf
done
