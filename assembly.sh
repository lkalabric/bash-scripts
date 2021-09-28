# Montagem DENV

# $1 Número da análise passado na linha de comando
if [ $# -eq 0 ]; then
	echo "Falta o método para mapeamento!"
	echo "Sintáxe: ./assembly.sh <METODO: 1,2,3,4>"
	exit 0
fi

# 1 Genome assembly using minimap2-miniasm pipeline (gera unitigs sequences)

if [ $1 -eq 1 ]; then
	# Fonte: https://timkahlke.github.io/LongRead_tutorials/ASS_M.html
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
	nanopolish index -d ../data/DENV_FTA_1/DENV_Run1_data/fast5_pass -s ../data/DENV_FTA_1/DENV_Run1_data/sequencing_summary/MT-110616_20190710_214507_FAK92171_minion_sequencing_run_DENV_FTA_1_sequencing_summary.txt albacore_output.fastq
	# Computa uma nova sequencia consenso
	minimap2 -t 12 -ax map-ont ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq -o "barcode01.$1.axligned.sam"
	samtools sort "barcode01.$1.axligned.sam" -o "barcode01.$1.axligned.sorted.bam" -T reads.tmp -
	samtools index "barcode01.$1.axligned.sorted.bam"
	python3 nanopolish_makerange.py ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta | parallel --results nanopolish.results -P 8 \
    nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq -b "barcode01.$1.axligned.sorted.bam" -g ../data/REFSEQ/Flaviviridae/NC_001477.1_DENV1.fasta -t 4 --min-candidate-frequency 0.1

#	nanopolish variantes --consensus  -d ../ngs-analysis/DENV_FTA_1_hac/wf31/PRINSEQ/barcode01.good.fastq -o "barcode01.$1.axligned.sam"
#	samtools sort "barcode01.$1.axligned.sam" -o "barcode01.$1.axligned.sorted.bam"
	exit 4
fi
