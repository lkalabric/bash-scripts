#!/bin/bash

# script: assembly2.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# última atualização: MAR, 18 2022
# versão 2: pipeline minimap2-samtools

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=$3		#Worflow de bioinformatica

if [[ $# -eq 0 ]]; then
	echo "Falta o nome da library,modelo Guppy Basecaller ou número do workflow!"
	echo "Sintáxe: ./assembly2.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,31>"
	exit 0
fi

# Caminho INPUT dos genomas referencia
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

# Parâmetro de otimização das análises por CPU
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros para controle da qualidade mínima
QSCORE=9
LENGTH=100

# Cria o arquivo índice das sequencias referencias para mapeamento das reads pelo minimap2
echo "Criando arquivos índices para os genomas referencia..."
for j in $(find ${REFGENDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
  REFGENFASTA="${REFGENDIR}/${j}.fasta"
  REFGENMMI="${REFGENDIR}/${j}.mmi"
  [ ! -f $REFGENMMI ] && minimap2 -d $REFGENMMI $REFGENFASTA
done  

# Mapeamento das sequencias em genomas referência e análise de cobertura
for j in $(find ${REFGENDIR} -type f -name "*.mmi" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		echo "Mapeando ${READSLEVELDIR}/${i}.corrected.fasta..."
		minimap2 -t ${THREADS} -ax map-ont ${REFGENDIR}/${j}.mmi ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.${j}.sorted.bam -	
		samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.${j}.sorted.bam > ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam
		samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam
		samtools mpileup -A -d 0 -Q 0 ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.${j}
	done
done

exit
