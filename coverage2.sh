#!/bin/bash

# script: coverage2.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# última atualização: MAR, 17 2022
# versão 2: roda samtools coverage & fastcov.py

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=$3		# Worflow de bioinformatica
TAXON=$4	# Análisa cobertura por taxon

if [[ $# -eq 0 ]]; then
	echo "Falta o nome da library,modelo Guppy Basecaller ou número do workflow!"
	echo "Sintáxe: ./coverage2.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,31>"
	exit 0
fi

# Caminhos de INPUT dos dados
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
READSLEVELDIR="${RESULTSDIR}/wf${WF}/READS_LEVEL"
ASSEMBLYDIR="${RESULTSDIR}/wf${WF}/ASSEMBLY"
REFGENDIR="${HOME}/data/REFGEN"

# Caminhos de OUTPUT das análises
COVERAGEDIR="${RESULTSDIR}/wf${WF}/COVERAGE"
[ ! -d $COVERAGEDIR ] && mkdir -vp $COVERAGEDIR

# Parâmetro de otimização das análises por CPU
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros para controle da qualidade mínima
QSCORE=9
LENGTH=100

# Cria o arquivo índice das sequencias referencias para mapeamento das reads pelo minimap2
echo "Gerando arquivos índices para os genomas referencia..."
for j in $(find ${REFGENDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
  REFGENFASTA="${REFGENDIR}/${j}.fasta"
  REFGENMMI="${REFGENDIR}/${j}.mmi"
  [ ! -f $REFGENMMI ] && minimap2 -d $REFGENMMI $REFGENFASTA
done  

# Análise de cobertura
source activate ngs
for j in $(find ${REFGENDIR} -type f -name "*.mmi" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		echo "Analisando a cobertura ${READSLEVELDIR}/${i}.corrected.fasta..."
		samtools coverage ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam > ${COVERAGEDIR}/${i}.${j}.coverage.txt
		samtools coverage -A -w 32 ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam > ${COVERAGEDIR}/${i}.${j}.histogram.txt
		samtools depth ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam > ${COVERAGEDIR}/${i}.${j}.depth.txt
		fastcov.py ${ASSEMBLYDIR}/${i}.${j}.sorted.mapped.bam -o ${COVERAGEDIR}/${i}.${j}.fastcov.pdf -c ${COVERAGEDIR}/${i}.${j}.fastcov.txt
	done
done

# Análise de cobertura da biblioteca por taxon
for j in $(find ${REFGENDIR} -type f -name "*.sorted.mapped.bam" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	echo "Determinando a cobertura da biblioteca para o taxon ${j}..."
	fastcov.py ${ASSEMBLYDIR}/barcode*.${j}.sorted.mapped.bam -o ${COVERAGEDIR}/${RUNNAME}_${MODEL}_${j}.fastcov.pdf
	fastcov.py -l ${ASSEMBLYDIR}/barcode*.${j}.sorted.mapped.bam -o ${COVERAGEDIR}/${RUNNAME}_${MODEL}_${j}.fastcov_log.pdf
done

exit
