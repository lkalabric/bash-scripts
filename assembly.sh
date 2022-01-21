#!/bin/bash

# script: assembly.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# objetivo: testar diferentes montadores
# criação: 25 AGO 2021
# ultima atualização: 17 OUT 2021
# atualização: configuração de variáveis e teste do método 6

OUTDIR=/home/brazil1/assembly
#REFSEQ="../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta"
REFSEQ="../data/REFSEQ/Togaviridae/NC_004162_CHIKV-S27.fasta"
#RUNNAME=DENV_FTA_1_hac
RUNNAME="NGS_LIBRARY13_fast"
BARCODE="barcode01"
SAMPLE="$HOME/ngs-analysis/$RUNNAME/wf31/PRINSEQ/$BARCODE.good.fastq"



# Valiação da entrada de dados na linha de comando
# $1 Número da análise passado na linha de comando
if [ $# -eq 0 ]; then
	echo "Falta o método para montagem!"
	echo "Sintáxe: ./assembly.sh <METODO: 1..6>"
	echo "Exemplo: ./assembly.sh 6"
	exit 0
fi


# 1 Genome assembly using minimap2-miniasm pipeline (gera unitigs sequences)

if [ $1 -eq 1 ]; then
	minimap2 -x ava-ont \
	 ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq \
	 ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq \
	| gzip -1 > "minimap.$1.paf.gz"

	miniasm -f \
	 ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq \
	" minimap.$1.paf.gz" > "miniasm.$1.gfa"
	awk '/^S/{print ">"$2"\n"$3}' "miniasm.$1.gfa" > "miniasm.$1.fasta"
	exit 1
fi

# 2 Mapea as reads usando um genoma referência
if [ $1 -eq 2 ]; then
	source activate ngs
	# long sequences against a reference genome
	minimap2 -t 12 -a ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq -o "barcode01.$1.mapped.sam"
	samtools sort "barcode01.$1.mapped.sam" -o "barcode01.$1.mapped.sorted.bam"
	fastcov -s "barcode01.$1.mapped.sorted.bam" -o "barcode01.$1.fastcov.pdf"
	exit 3
fi


# 3 Mapea as reads usando um genoma referência
if [ $1 -eq 3 ]; then
	source activate ngs
	# Cria um indice antes de mapear
	minimap2 -t 12 -a ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq -o "barcode01.$1.mapped.sam"
	samtools sort "barcode01.$1.mapped.sam" -o "barcode01.$1.mapped.sorted.bam"
	exit 3
fi

# 4 Mapea as reads usando um genoma referência
if [ $1 -eq 4 ]; then
	source activate ngs
	# use presets (no test data) # Oxford Nanopore genomic reads
	minimap2 -t 12 -ax map-ont ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq -o "barcode01.$1.axligned.sam"
	samtools sort "barcode01.$1.axligned.sam" -o "barcode01.$1.axligned.sorted.bam"
	exit 4
fi

# 5 Montagem da sequencia consenso usando um genoma referência
if [ $1 -eq 5 ]; then
	# Fonte: https://github.com/jts/nanopolish
	# Pré-processamento dos dados
	nanopolish index -d ../data/DENV_FTA_1/DENV_Run1_data/fast5_pass/ -s ../data/DENV_FTA_1/DENV_Run1_data/sequencing_summary/MT-110616_20190710_214507_FAK92171_minion_sequencing_run_DENV_FTA_1_sequencing_summary.txt "${OUTDIR}/barcode01.fasta"
	# Computa uma nova sequencia consenso
	minimap2 -t 12 -ax map-ont ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq -o "$OUTDIR/barcode01.$1.axligned.sam"
	samtools sort "$OUTDIR/barcode01.$1.axligned.sam" -o "$OUTDIR/barcode01.$1.axligned.sorted.bam"
	samtools index "$OUTDIR/barcode01.$1.axligned.sorted.bam"
	# Quebra o genoma em pedações de 50Kb e monta em paralelo
	#	python3 nanopolish_makerange.py ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta | parallel --results nanopolish.results -P 8 \
	nanopolish variants -o "${OUTDIR}/polished.vcf" -r "$OUTDIR/barcode01.fasta" -b "$OUTDIR/barcode01.$1.axligned.sorted.bam" -g ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta -t 4 --min-candidate-frequency 0.1 -p 1
#	nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq -b "barcode01.$1.axligned.sorted.bam" -g ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta -t 4 --min-candidate-frequency 0.1

#	nanopolish variantes --consensus  -d ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq -o "barcode01.$1.axligned.sam"
#	samtools sort "barcode01.$1.axligned.sam" -o "barcode01.$1.axligned.sorted.bam"
	exit 5
fi

# 6 Montagem da sequencia consenso usando um genoma referência
if [ $1 -eq 6 ]; then
	# Fonte: https://github.com/jts/nanopolish
	# Indexando a sequencia referencia
	bwa index $REFSEQ
	bwa mem $REFSEQ "../ngs-analysis/$RUNNAME/wf31/PRINSEQ/$BARCODE.good.fastq" > "$OUTDIR/$BARCODE.$1.bwa-mem.sam"
	samtools sort "$OUTDIR/$BARCODE.$1.bwa-mem.sam" -o "$OUTDIR/$BARCODE.$1.bwa-mem.sorted.bam"
	samtools index "$OUTDIR/$BARCODE.$1.bwa-mem.sorted.bam"
	samtools coverage "$OUTDIR/$BARCODE.$1.bwa-mem.sorted.bam" -m -o "$OUTDIR/$BARCODE.$1.coverage"
	cat "../assembly/$BARCODE.6.coverage"
	fastcov.py "$OUTDIR/$BARCODE.$1.bwa-mem.sorted.bam" -o "../assembly/$BARCODE.6.fastcov.pdf"
	exit 6
fi

# 7 Montagem utilizando wtdbg2
if [ $1 -eq 7 ]; then
	# Fonte: https://github.com/ruanjue/wtdbg2
	
	# Ativar o ambiente Conda
	soruce activate ngs
	
	# assemble long reads
	wtdbg2 -x ont -g 4.6m -i $SAMPLE -t 12 -fo dbg

	# derive consensus
	wtpoa-cns -t 16 -i dbg.ctg.lay.gz -fo dbg.raw.fa

	# polish consensus, not necessary if you want to polish the assemblies using other tools
	#minimap2 -t16 -ax map-pb -r2k dbg.raw.fa reads.fa.gz | samtools sort -@4 >dbg.bam
	#samtools view -F0x900 dbg.bam | ./wtpoa-cns -t 16 -d dbg.raw.fa -i - -fo dbg.cns.fa

	# Addtional polishment using short reads
	#bwa index dbg.cns.fa
	#bwa mem -t 16 dbg.cns.fa sr.1.fa sr.2.fa | samtools sort -O SAM | ./wtpoa-cns -t 16 -x sam-sr -d dbg.cns.fa -i - -fo dbg.srp.fa
	exit 7
fi
