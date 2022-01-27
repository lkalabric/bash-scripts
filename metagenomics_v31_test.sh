#!/bin/bash

# script: metagenomics_v31_test.sh
# autores: Luciano Kalabric <luciano.kalabric@fiocruz.br> & Alessandra Gonzalez
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 27 JAN 2022
# última atualização: 27 JAN 2022

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=$3		# Workflow de bioinformatica de teste

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
DEMUXDIR="${RESULTSDIR}/test${WF}/DEMUX"
DEMUXCATDIR="${RESULTSDIR}/test${WF}/DEMUX_CAT"
CUTADAPTDIR="CUTADAPT"
NANOFILTDIR="${RESULTSDIR}/test${WF}/NANOFILT"
PRINSEQDIR="${RESULTSDIR}/test${WF}/PRINSEQ"
QUERYDIR="${RESULTSDIR}/test${WF}/QUERY"
BLASTDIR="${RESULTSDIR}/test${WF}/BLAST"
READSLEVELDIR="${RESULTSDIR}/test${WF}/READS_LEVEL"
CONTIGLEVELDIR="${RESULTSDIR}/test${WF}/CONTIGS_LEVEL"
ASSEMBLYDIR="${RESULTSDIR}/test${WF}/ASSEMBLY"

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

# Sumario da corrida (dados disponíveis no arquivo report*.pdf)
echo "Sumário da corrida"
echo "Total files:"
ls $(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;) | wc -l
echo "Total reads:"
# h5ls "$(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;)"/*.fast5 | wc -l

# Step 1 - Basecalling
# Esta etapa pode ser realizada pelo script guppy_gpu_v1_ag.sh no LAPTOP-Yale
if [ ! -d $BASECALLDIR ]; then
	echo -e "\nExecutando guppy_basecaller..."
	mkdir -vp $BASECALLDIR
	# Step 1 - Basecalling (comum a todos workflows)
	# Esta etapa está sendo realizada pelo script guppy_gpu_v1_ag.sh no LAPTOP-Yale
	# Comando para guppy_basecaller usando GPU
	guppy_basecaller -r -i ${RAWDIR} -s "${BASECALLDIR}" -c ${CONFIG} -x auto  --gpu_runners_per_device ${GPUPERDEVICE} --chunk_size ${CHUNCKSIZE} --chunks_per_runner ${CHUNKPERRUNNER} --verbose_logs
fi

#
# Test 1
#
if [[ $WF -eq 1 ]]; then 
	# Step 2 - Demultiplex & adapter removal
	echo -e "\nExecutando guppy_barcoder.1..."
	# Sem headcrop 18 (uso cutadapt)
	mkdir -vp $DEMUXDIR.1
	guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --trim_barcodes  

	# Move a pasta contendo as reads unclassified para barcode00
	[ -d "${DEMUXDIR}/unclassified" ] && mv "${DEMUXDIR}/unclassified" "${DEMUXDIR}/barcode00"

	# Concatena todos arquivos .fastq de cada barcode em um arquivo .fastq único
	[ ! -d ${DEMUXCATDIR} ] && mkdir -vp ${DEMUXCATDIR}
	for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
	    [ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > "${DEMUXCATDIR}/${i}.fastq"
	done

	# Step 2.1 - Remoção dos primers
	echo -e "\nExecutando cutadapt..."
	[ ! -d ${DEMUXCATDIR}/${CUTADAPTDIR} ] && mkdir -vp ${DEMUXCATDIR}/${CUTADAPTDIR}
	for i in $(find ${DEMUXCATDIR} -type f -exec basename {} .fastq \;); do
		cutadapt -g ${PRIMER} -e 0.2 --discard-untrimmed -o "${DEMUXCATDIR}/${CUTADAPTDIR}/${i}.fastq" "${DEMUXCATDIR}/${i}.fastq"
		echo -e "\nResultados ${i} $(grep -c "runid" ${DEMUXCATDIR}/${CUTADAPTDIR}/${i}.fastq | cut -d : -f 2 | awk '{s+=$1} END {printf "%.0f\n",s}')"
	done
	mv ${DEMUXCATDIR}/${CUTADAPTDIR}/*.fastq ${DEMUXCATDIR}
fi

#
# Test 2
#
if [[ $WF -eq 2 ]]; then 
	# Step 2 - Demultiplex & adapter removal
	echo -e "\nExecutando guppy_barcoder.2..."
	# Com headcrop 18
	mkdir -vp $DEMUXDIR
	guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER}

	# Move a pasta contendo as reads unclassified para barcode00
	[ -d "${DEMUXDIR}/unclassified" ] && mv "${DEMUXDIR}/unclassified" "${DEMUXDIR}.2/barcode00"

	[ ! -d ${DEMUXCATDIR}.2 ] && mkdir -vp ${DEMUXCATDIR}.2
	for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
	    [ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > "${DEMUXCATDIR}/${i}.fastq"
	done
fi

# WF 3 - Classificação Taxonômica pelo Kraken2

	source activate ngs
	# Step 4 - Filtro por tamanho
	echo -e "\nExecutando NanoFilt..."
	[ ! -d ${NANOFILTDIR} ] && mkdir -vp ${NANOFILTDIR}
	for i in $(find ${DEMUXCATDIR} -type f -exec basename {} .fastq \;); do
		NanoFilt -l ${LENGTH} < "${DEMUXCATDIR}/${i}.fastq" > "${NANOFILTDIR}/${i}.fastq" 
		# Resultados disponíveis no report do Prinseq (Input sequences) 
	done

	# Step 5 - Filtro de complexidade
	# Link: https://chipster.csc.fi/manual/prinseq-complexity-filter.html
	echo -e "\nExecutando prinseq-lite.pl..."
	[ ! -d ${PRINSEQDIR} ] && mkdir -vp ${PRINSEQDIR}
	for i in $(find ${NANOFILTDIR} -type f -exec basename {} .fastq \;); do
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

echo "Workflow $WF concluido com sucesso!"
exit 0

