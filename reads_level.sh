#!/bin/bash

# script: reads_level.sh
# autores: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 14 ABR 2022
# última atualização: 14 ABR 2022

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=$3		# Workflow de bioinformatica 1, 2 ou 31

if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, número do worflow ou modelo Guppy Basecaller!"
	echo "Sintáxe: ./reads_level.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,31>"
	exit 0
fi

# Caminho de INPUT dos bancos de dados
REFSEQDIR=${HOME}/data/REFSEQ
HUMANREFDIR=${HOME}/data/GRCh38
BLASTDBDIR="${HOME}/data/BLAST_DB"
KRAKENDB="${HOME}/data/KRAKEN2_DB" # Substituir pelo nosso banco de dados se necessário KRAKEN2_USER_DB
echo "Preparando pastas para (re-)análise dos dados..."
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
DEMUXCATDIR="${RESULTSDIR}/DEMUX_CAT"
FILTER_BY_START_TIMEDIR="${RESULTSDIR}/FILTER_BY_START_TIME"

# Caminhos de OUTPUT das análises
CUTADAPTDIR="${RESULTSDIR}/wf${WF}_filter/CUTADAPT"
NANOFILTDIR="${RESULTSDIR}/wf${WF}_filter/NANOFILT"
PRINSEQDIR="${RESULTSDIR}/wf${WF}_filter/PRINSEQ"
QUERYDIR="${RESULTSDIR}/wf${WF}_filter/QUERY"
BLASTDIR="${RESULTSDIR}/wf${WF}_filter/BLAST"
READSLEVELDIR="${RESULTSDIR}/wf${WF}_filter/READS_LEVEL"
CONTIGLEVELDIR="${RESULTSDIR}/wf${WF}_filter/CONTIGS_LEVEL"
ASSEMBLYDIR="${RESULTSDIR}/wf${WF}_filter/ASSEMBLY"

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

# WF 1 - Classificação Taxonômica pelo Epi2ME
# Copiar a pasta /pass para o Epi2ME
if [[ $WF -eq 1 ]]; then
	exit 1
fi

# if false; then # Desvio para execução rápida

# WF 2 & 31
#
# Read-level taxid
#

# WF 2 - Classificação Taxonômica através de busca no BLASTDB local
if [[ $WF -eq 2 ]]; then 
	# Step 4 - Remoção dos primers
	echo -e "\nExecutando cutadapt..."
	[ ! -d ${CUTADAPTDIR} ] && mkdir -vp ${CUTADAPTDIR}
	for i in $(find ${FILTER_BY_START_TIMEDIR} -type f -exec basename {} .fastq \; | sort); do
		cutadapt -g ${PRIMER} -e 0.2 --discard-untrimmed -o "${CUTADAPTDIR}/${i}.fastq" "${DEMUXCATDIR}/${i}.fastq"
		echo -e "\nResultados ${i} $(grep -c "runid" ${CUTADAPTDIR}/${i}.fastq | cut -d : -f 2 | awk '{s+=$1} END {printf "%.0f\n",s}')"
	done

	# Step 5 - Filtro por tamanho
	echo -e "\nExecutando NanoFilt..."
	[ ! -d ${NANOFILTDIR} ] && mkdir -vp ${NANOFILTDIR}
	for i in $(find "${CUTADAPTDIR}" -type f -exec basename {} .fastq \; | sort); do
		NanoFilt -l ${LENGTH} < "${CUTADAPTDIR}/${i}.fastq" > "${NANOFILTDIR}/${i}.fastq" 
		# Resultados disponíveis no report do Prinseq (Input sequences) 
	done

	# Step 6 - Filtro de complexidade
	# Link: https://chipster.csc.fi/manual/prinseq-complexity-filter.html
	echo -e "\nExecutando prinseq-lite.pl..."
	[ ! -d ${PRINSEQDIR} ] && mkdir -vp ${PRINSEQDIR}
	for i in $(find ${NANOFILTDIR} -type f -exec basename {} .fastq \; | sort); do
		echo -e "\nResultados ${i}..."
		prinseq-lite.pl -fastq "${NANOFILTDIR}/${i}.fastq" -out_good "${PRINSEQDIR}/${i}.good" -out_bad "${PRINSEQDIR}/${i}.bad -graph_data" "${PRINSEQDIR}/${i}.gd" -no_qual_header -lc_method dust -lc_threshold 40
		# Resultados disponíveis no report do Prinseq (Good sequences)
	done

	# Converte arquivos .fastq em .fasta para query no blastn
	[ ! -d "${QUERYDIR}" ] && mkdir -vp ${QUERYDIR}
	for i in $(find "${PRINSEQDIR}"/*.good.fastq -type f -exec basename {} .good.fastq \; | sort); do
		sed -n '1~4s/^@/>/p;2~4p' "${PRINSEQDIR}/${i}.good.fastq" > "${QUERYDIR}/${i}.fasta"
	done
		
	# Step 7 - Classificação taxonômica utilizando blastn
	# Preparação do BLASTDB local
	# Script: makeblastdb_refseq.sh
		# Concatena todas as REFSEQs num arquivo refseq.fasta único e cria o BLASTDB
		# Extrai do arquvio refseq.fasta a lista acesso refseq.acc
		# Cria a partir do arquivo refseq.acc o arquivo refseq.map que mapeia os taxid (números que identificam as espécies taxonômica)
	# Busca as QUERIES no BLASTDB local
	echo -e "\nExecutando blastn..."
	[ ! -d ${BLASTDIR} ] && mkdir -vp ${BLASTDIR}
	for i in $(find ${QUERYDIR} -type f -exec basename {} .fasta \; | sort); do
		echo -e "\nAnalisando dados ${BLASTDIR}/${i}..."
		blastn -db "${BLASTDBDIR}/refseq" -query "${QUERYDIR}/${i}.fasta" -out "${BLASTDIR}/${i}.blastn" -outfmt "6 sacc staxid" -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
		# Busca remota
		# blastn -db nt -remote -query ${QUERYDIR}/${i}.fasta -out ${BLASTDIR}/${i}.blastn -outfmt "6 qacc saccver pident sscinames length mismatch gapopen evalue bitscore"  -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
		echo -e "\nResultados ${i}"
		~/scripts/blast_report.sh "${RUNNAME}_${MODEL}" "${i}"
	done
	exit 2
fi

# WF 3 - Classificação Taxonômica pelo Kraken2
if [[ $WF -eq 31 ]]; then 
	source activate ngs
	# Step 4 - Filtro por tamanho
	echo -e "\nExecutando NanoFilt..."
	[ ! -d ${NANOFILTDIR} ] && mkdir -vp ${NANOFILTDIR}
	for i in $(find ${FILTER_BY_START_TIMEDIR} -type f -exec basename {} .fastq \; | sort); do
		NanoFilt -l ${LENGTH} < "${FILTER_BY_START_TIMEDIR}/${i}.fastq" > "${NANOFILTDIR}/${i}.fastq" 
		# Resultados disponíveis no report do Prinseq (Input sequences) 
	done

	# Step 5 - Filtro de complexidade
	# Link: https://chipster.csc.fi/manual/prinseq-complexity-filter.html
	echo -e "\nExecutando prinseq-lite.pl..."
	[ ! -d ${PRINSEQDIR} ] && mkdir -vp ${PRINSEQDIR}
	for i in $(find ${NANOFILTDIR} -type f -exec basename {} .fastq \; | sort); do
		echo -e "\nCarregando os dados ${i}..."
		prinseq-lite.pl -fastq "${NANOFILTDIR}/${i}.fastq" -out_good "${PRINSEQDIR}/${i}.good" -out_bad "${PRINSEQDIR}/${i}.bad -graph_data" "${PRINSEQDIR}/${i}.gd" -no_qual_header -lc_method dust -lc_threshold 40
		# Resultados disponíveis no report do Prinseq (Good sequences)
	done

	# Step 6 - Remoção das reads do genoma humano
	echo -e "\nExecutando minimap2 & samtools para filtrar as reads do genoma humano..."
	[ ! -d "${READSLEVELDIR}" ] && mkdir -vp ${READSLEVELDIR}
	[ ! -d "${ASSEMBLYDIR}" ] && mkdir -vp ${ASSEMBLYDIR}
	# Cria o arquivo índice do genoma humano para reduzir o tempo de alinhamento
	minimap2 -d ${HUMANREFMMI} ${HUMANREFSEQ}
	# Loop para analisar todos barcodes, um de cada vez
	for i in $(find ${PRINSEQDIR} -type f -name "*.good.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		echo -e "\nCarregando os dados ${i}..."
	    	# Alinha as reads contra o arquivo indice do genoma humano e ordena os segmentos
	    	minimap2 -ax map-ont -t ${THREADS} ${HUMANREFMMI} ${PRINSEQDIR}/${i}.good.fastq | samtools sort -@ 12 -o ${READSLEVELDIR}/${i}.sorted.bam -
	    	# Indexa o arquivo para acesso mais rápido
	    	samtools index -@ ${THREADS} ${READSLEVELDIR}/${i}.sorted.bam
	    	# Filtra os segmento não mapeados Flag 4 (-f 4) para um novo arquivo .bam 
	    	samtools view -bS -f 4 ${READSLEVELDIR}/${i}.sorted.bam > ${READSLEVELDIR}/${i}.unmapped.sam -@ 12
		# Salva os dados no formato .fastq
		samtools fastq -f 4 ${READSLEVELDIR}/${i}.unmapped.sam > ${READSLEVELDIR}/${i}.unmapped.fastq -@ 12
	done

	# Step 7 - Autocorreção das reads
	echo -e "\nExecutando minimap2 & racon para autocorreção das reads contra a sequencia consenso..."
	for i in $(find ${PRINSEQDIR} -type f -name "*.good.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		echo -e "\nCarregando os dados ${i}..."
		# Alinhar todas as reads com elas mesmas para produzir sequencias consenso a partir do overlap de reads
		minimap2 -ax ava-ont -t ${THREADS} ${READSLEVELDIR}/${i}.unmapped.fastq ${READSLEVELDIR}/${i}.unmapped.fastq > ${READSLEVELDIR}/${i}.overlap.sam
		# Correção de erros a partir das sequencias consenso
	 	racon -t ${THREADS} -f -u ${READSLEVELDIR}/${i}.unmapped.fastq ${READSLEVELDIR}/${i}.overlap.sam ${READSLEVELDIR}/${i}.unmapped.fastq > ${READSLEVELDIR}/${i}.corrected.fasta
	done

	# Step 8 - Classificação taxonômica
	echo -e "\nExecutando o Kraken2..."
	for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		# kraken2 --db ${KRAKENDB} --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --report-minimizer-data --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta
		echo -e "\nCarregando os dados ${i}..."
		kraken2 --db ${KRAKENDB} --quick --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta
		echo -e "\nResultados ${i}"
		~/scripts/kraken2_quick_report.sh "${RUNNAME}_${MODEL}" "${i}"
	done
	exit 31
fi

echo "Workflow $WF concluido com sucesso!"
exit 0
