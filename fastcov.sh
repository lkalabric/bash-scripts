#!/bin/bash

# script: metagenomics_v4_lk.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# última atualização: SET, 17 2021
# versão 4: roda guppy_basecaller e guppy_barcoder apenas uma vez

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup

if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, número do worflow ou modelo Guppy Basecaller!"
	echo "Sintáxe: ./metagenomic_v3_lk.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,3>"
	exit 0
fi

# Caminho INPUT dos bancos de dados
REFSEQDIR="${HOME}/data/REFSEQ"
HUMANREFDIR="${HOME}/data/GRCh38"
BLASTDBDIR="${HOME}/data/BLAST_DB"
KRAKENDB="${HOME}/data/KRAKEN2_DB" # Substituir pelo nosso banco de dados se necessário KRAKEN2_USER_DB
TESTEDIR="${HOME}/data/teste"

# Caminhos de OUTPUT das análises
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
BASECALLDIR="${RESULTSDIR}/BASECALL"
DEMUXDIR="${RESULTSDIR}/DEMUX"
DEMUXCATDIR="${RESULTSDIR}/DEMUX_CAT"
CUTADAPTDIR="${RESULTSDIR}/CUTADAPT"
NANOFILTDIR="${RESULTSDIR}/NANOFILT"
PRINSEQDIR="${RESULTSDIR}/PRINSEQ"
QUERYDIR="${RESULTSDIR}/QUERY"
BLASTDIR="${RESULTSDIR}/BLAST"
READSLEVELDIR="${RESULTSDIR}/READS_LEVEL"
CONTIGSLEVELDIR="${RESULTSDIR}/CONTIGS_LEVEL"
ASSEMBLYDIR="${RESULTSDIR}/ASSEMBLY"
FASTCOVDIR="${RESULTSDIR}/FASTCOV"

# Parâmetro de otimização das análises por CPU
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros para controle da qualidade mínima
QSCORE=9
LENGTH=100

# Parâmetros minimap2
# Mapeamento contra genomas referência e plot de cobertura
CHIKVREFSEQ="${TESTEDIR}/NC_004162_CHIKV-S27.fasta"
CHIKVREFMMI="${TESTEDIR}/NC_004162_CHIKV-S27.mmi"

[ ! -f "${TESTEDIR}/NC_004162.2_CHIKV-S27.mmi" ] && minimap2 -d $CHIKVREFMMI $CHIKVREFSEQ



DENV1REFSEQ="${REFSEQDIR}/Flaviviridae/NC_001477.1_DENV1.fasta"
DENV2REFSEQ="${REFSEQDIR}/Flaviviridae/NC_001474.2_DENV2.fasta"
DENV3REFSEQ="${REFSEQDIR}/Flaviviridae/NC_001475.2_DENV3.fasta"
DENV4REFSEQ="${REFSEQDIR}/Flaviviridae/NC_002640.1_DENV4.fasta"
ZIKVREFSEQ="${REFSEQDIR}/Flaviviridae/NC_012532.1_ZIKV.fasta"

[ ! -d $ASSEMBLYDIR ] && mkdir -vp $ASSEMBLYDIR
[ ! -d $FASTCOVDIR ] && mkdir -vp $FASTCOVDIR


# Mapeamento CHIKV
for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	echo "Mapeando ${READSLEVELDIR}/${i}.corrected.fasta..."
	minimap2 -t ${THREADS} -ax map-ont ${CHIKVREFMMI} ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.chikv.sorted.bam -	
	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.chikv.sorted.bam > ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
	samtools mpileup -A -d 0 -Q 0 ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.chikv
	
done

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
	echo "Determinando a cobertura do assembly ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam..."
	fastcov ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam -o ${FASTCOVDIR}/${i}.assembly_chikv.pdf
	fastcov -l -o ${FASTCOVDIR}/${i}.assembly_chikv_log.png {ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam 
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
done
