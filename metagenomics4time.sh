#!/bin/bash

# script: metagenomics4time.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 09 JUN 2022
# última atualização: 09 JUN 2022
# versão 1: modulariza as etapas do workflow e permite criar diferentes wokflows executado cada etapa como uma função

# Descrição de cada etapa disponível para construção dos workflows
# sequencing_summary1
# sequencing_summary2
# basecalling
# demux_cat1
# demux_cat2
# qc_filter1
# qc_filter2
# human_filter
# autocorrection
# blast
# kraken

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=$3		# Workflow de bioinformatica 1, 2 ou 3
if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, número do worflow ou modelo Guppy Basecaller!"
	echo "Sintáxe: ./metagenomics4.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,3>"
	exit 0
fi

# Caminho de INPUT dos dados fast5
RAWDIR="${HOME}/data/${RUNNAME}" # Se análise começar com o Barcoder
if [ ! -d $RAWDIR ]; then
	echo "Pasta de dados não encontrada!"
	exit 0
fi

# Caminho de INPUT dos bancos de dados
HUMANREFDIR=${HOME}/data/GRCh38
REFSEQDIR=${HOME}/data/REFSEQ
BLASTDBDIR="${HOME}/data/BLAST_DB"
KRAKENDB="${HOME}/data/KRAKEN2_DB" # Substituir pelo nosso banco de dados se necessário KRAKEN2_USER_DB

# Caminhos de OUTPUT das análises
echo "Preparando pastas para (re-)análise dos dados..."
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}/time"
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

# Parâmetro de otimização minimap2, samtools, racon e kraken2
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros de qualidade mínima
  QSCORE=9
  LENGTH=100

# Parâmetros Barcoder ou Cutadapt para remoção do primer
TRIMADAPTER=18
PRIMER="GTTTCCCACTGGAGGATA"

# Pausa a execução para debug
# read -p "Press [Enter] key to continue..."

function sequencing_summary1 () {
  $RAWDIR=$1
  # Sumario do sequenciamento (dados disponíveis no arquivo report*.pdf)
  echo "Sumário da corrida"
  echo "Total files:"
  ls $(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;) | wc -l
  echo "Total reads:"
  # h5ls "$(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;)"/*.fast5 | wc -l
}

function basecalling () {
  # Basecalling  (comum a todos workflows)
  $RAWDIR=$1
  $BASECALLDIR=$2
  $MODEL=$3
  # Parâmetros Guppy basecaller (ONT)
  CONFIG="dna_r9.4.1_450bps_${MODEL}.cfg" #dna_r9.4.1_450bps_fast.cfg dna_r9.4.1_450bps_hac.cfg dna_r9.4.1_450bps_sup.cfg
  # Parâmetros para otimização do Guppy basecaller para modelo fast utilizadando o LAPTOP-Yale (benckmark)
  case $MODEL in
    fast)
      GPUPERDEVICE=4		
      CHUNCKSIZE=1000		
      CHUNKPERRUNNER=50
    ;;
    hac)
      GPUPERDEVICE=12		
      CHUNCKSIZE=2000		
      CHUNKPERRUNNER=256
    ;;
    sup)
      GPUPERDEVICE=12		
      CHUNCKSIZE=1000		
      CHUNKPERRUNNER=256
    ;;
    *)
      GPUPERDEVICE=4		
      CHUNCKSIZE=1000		
      CHUNKPERRUNNER=50
    ;;
  esac
  echo "Parâmetros otimizados para guppy_basecaller no modelo selecionado: $MODEL"
  echo "GPUPERDEVICE: $GPUPERDEVICE"
  echo "CHUNCKSIZE: $CHUNCKSIZE"
  echo "CHUNKPERRUNNER: $CHUNKPERRUNNER"
  if [ ! -d $BASECALLDIR ]; then
	echo -e "\nExecutando guppy_basecaller..."
	mkdir -vp $BASECALLDIR
	# Esta etapa está sendo realizada pelo script guppy_gpu_v1_ag.sh no LAPTOP-Yale
	# Comando para guppy_basecaller usando GPU
	guppy_basecaller -r -i ${RAWDIR} -s "${BASECALLDIR}" -c ${CONFIG} -x auto  --gpu_runners_per_device ${GPUPERDEVICE} --chunk_size ${CHUNCKSIZE} --chunks_per_runner ${CHUNKPERRUNNER} --verbose_logs
  fi
}

# Para debug
# if false; then # Desvio para execução rápida
# fi # Fim do desvio para execução rápida

function demux_cat1 () {
  # Demultiplex, adapter removal & sem headcrop 18
  $BASECALLDIR=$1
  $DEMUXDIR=$2
  $DEMUXCATDIR=$3
  # Parâmetros Guppy barcoder (ONT)
  ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
  if [ ! -d $DEMUXDIR ]; then
    echo -e "\nExecutando guppy_barcoder..."
    mkdir -vp $DEMUXDIR
    #bm5 WF2 sem headcrop 18 (uso cutadapt)
    guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --trim_barcodes  
  fi
  # Move a pasta contendo as reads unclassified para barcode00
  [ -d "${DEMUXDIR}/unclassified" ] && mv "${DEMUXDIR}/unclassified" "${DEMUXDIR}/barcode00"

  # Concatena todos arquivos .fastq de cada barcode em um arquivo .fastq único
  [ ! -d ${DEMUXCATDIR} ] && mkdir -vp ${DEMUXCATDIR}
  for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
      [ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > "${DEMUXCATDIR}/${i}.fastq"
  done
}

function demux_cat2 () {
  # Demultiplex, adapter removal com headcrop 18
  $BASECALLDIR=$1
  $DEMUXDIR=$2
  $DEMUXCATDIR=$3
  # Parâmetros Guppy barcoder (ONT)
  ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
  if [ ! -d $DEMUXDIR ]; then
    echo -e "\nExecutando guppy_barcoder..."
    mkdir -vp $DEMUXDIR
    #bm7 WF31 com headcrop 18
    guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER}
  fi
  # Move a pasta contendo as reads unclassified para barcode00
  [ -d "${DEMUXDIR}/unclassified" ] && mv "${DEMUXDIR}/unclassified" "${DEMUXDIR}/barcode00"
  # Concatena todos arquivos .fastq de cada barcode em um arquivo .fastq único
  [ ! -d ${DEMUXCATDIR} ] && mkdir -vp ${DEMUXCATDIR}
  for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
      [ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > "${DEMUXCATDIR}/${i}.fastq"
  done
}

function sequencing_summary2 () {
	# pycoQC summary
	$RESULTSDIR=$1
	$RUNNAME=$2
	$BASECALLDIR=$3
	$DEMUXDIR=$4
	$QSCORE=$5
	$LENGTH=$6
	echo -e "\nExecutando pycoQC..."
	# source activate ngs
	# Default do guppybasecaller para min_pass_qual 8 (isso não está escrito em nenhum lugar)
	if [ ! -f "${RESULTSDIR}/${RUNNAME}_basecaller_pycoqc.html" ]; then
	    # Comando para pycoQC version 2.5
	    pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -o "${RESULTSDIR}/${RUNNAME}_basecaller_pycoqc.html" --report_title $RUNNAME --min_pass_qual 8
	fi
	if [ ! -f "${RESULTSDIR}/${RUNNAME}_pycoqc.html" ]; then
	    # Comando para pycoQC version 2.5
	    pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/barcoding_summary.txt" -o "${RESULTSDIR}/${RUNNAME}_pycoqc.html" --report_title $RUNNAME --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}
	fi
}

function primer_removal () {
	# Remoção dos primers
	$CUTADAPTDIR=$1
	$PRIMER=$2
	DEMUXCATDIR=$3
	echo -e "\nExecutando cutadapt..."
	[ ! -d ${CUTADAPTDIR} ] && mkdir -vp ${CUTADAPTDIR}
	for i in $(find ${DEMUXCATDIR} -type f -exec basename {} .fastq \;); do
		cutadapt -g ${PRIMER} -e 0.2 --discard-untrimmed -o "${CUTADAPTDIR}/${i}.fastq" "${DEMUXCATDIR}/${i}.fastq"
		echo -e "\nResultados ${i} $(grep -c "runid" ${CUTADAPTDIR}/${i}.fastq | cut -d : -f 2 | awk '{s+=$1} END {printf "%.0f\n",s}')"
	done
}

function qc_filter1 () {
	# Filtro por tamanho
	$CUTADAPTDIR=$1
	$LENGTH=$2
	$NANOFILTERDIR=$3
	$PRINSEQDIR=$4
	$QUERYDIR=$5
	echo -e "\nExecutando NanoFilt..."
	[ ! -d ${NANOFILTDIR} ] && mkdir -vp ${NANOFILTDIR}
	for i in $(find "${CUTADAPTDIR}" -type f -exec basename {} .fastq \;); do
		NanoFilt -l ${LENGTH} < "${CUTADAPTDIR}/${i}.fastq" > "${NANOFILTDIR}/${i}.fastq" 
		# Resultados disponíveis no report do Prinseq (Input sequences) 
	done

	# Filtro de complexidade
	# Link: https://chipster.csc.fi/manual/prinseq-complexity-filter.html
	echo -e "\nExecutando prinseq-lite.pl..."
	[ ! -d ${PRINSEQDIR} ] && mkdir -vp ${PRINSEQDIR}
	for i in $(find ${NANOFILTDIR} -type f -exec basename {} .fastq \;); do
		echo -e "\nResultados ${i}..."
		prinseq-lite.pl -fastq "${NANOFILTDIR}/${i}.fastq" -out_good "${PRINSEQDIR}/${i}.good" -out_bad "${PRINSEQDIR}/${i}.bad -graph_data" "${PRINSEQDIR}/${i}.gd" -no_qual_header -lc_method dust -lc_threshold 40
		# Resultados disponíveis no report do Prinseq (Good sequences)
	done

	# Converte arquivos .fastq em .fasta para query no blastn
	[ ! -d "${QUERYDIR}" ] && mkdir -vp ${QUERYDIR}
	for i in $(find "${PRINSEQDIR}"/*.good.fastq -type f -exec basename {} .good.fastq \;); do
		sed -n '1~4s/^@/>/p;2~4p' "${PRINSEQDIR}/${i}.good.fastq" > "${QUERYDIR}/${i}.fasta"
	done
}

function blast () {
	# Classificação taxonômica utilizando blastn
	$BLASTDIR=$1
	$QUERYDIR=$2
	$RUNNAME=$3
	$MODEL=$4
	# Preparação do BLASTDB local
	# Script: makeblastdb_refseq.sh
		# Concatena todas as REFSEQs num arquivo refseq.fasta único e cria o BLASTDB
		# Extrai do arquvio refseq.fasta a lista acesso refseq.acc
		# Cria a partir do arquivo refseq.acc o arquivo refseq.map que mapeia os taxid (números que identificam as espécies taxonômica)
	# Busca as QUERIES no BLASTDB local
	echo -e "\nExecutando blastn..."
	[ ! -d ${BLASTDIR} ] && mkdir -vp ${BLASTDIR}
	for i in $(find ${QUERYDIR} -type f -exec basename {} .fasta \;); do
		echo -e "\nAnalisando dados ${BLASTDIR}/${i}..."
		blastn -db "${BLASTDBDIR}/refseq" -query "${QUERYDIR}/${i}.fasta" -out "${BLASTDIR}/${i}.blastn" -outfmt "6 sacc staxid" -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
		# Busca remota
		# blastn -db nt -remote -query ${QUERYDIR}/${i}.fasta -out ${BLASTDIR}/${i}.blastn -outfmt "6 qacc saccver pident sscinames length mismatch gapopen evalue bitscore"  -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
		echo -e "\nResultados ${i}"
		~/scripts/blast_report.sh "${RUNNAME}_${MODEL}" "${i}"
	done
}

function qc_filter2 () {
	# Filtro por tamanho
	$DEMUXCATDIR=$1
	$LENGTH=$2
	$NANOFILTERDIR=$3
	$PRINSEQDIR=$4
	source activate ngs
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
}

function human_filter () {
	# Remoção das reads do genoma humano
	$HUMANREFDIR=$1
	$READSLEVELDIR=$2
	$ASSEMBLYDIR=$3
	$PRINSEQDIR=$4
	# Parâmetros minimap2 
	# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -P ${HUMANREFDIR}
	HUMANREFSEQ="${HUMANREFDIR}/GRCh38.p13.genome.fa.gz"
	HUMANREFMMI="${HUMANREFDIR}/GRCh38.p13.genome.mmi"
	# Cria o arquivo índice do genoma humano para reduzir o tempo de alinhamento
	# minimap2 -d ${HUMANREFMMI} ${HUMANREFSEQ}
	# Cria o arquivo índice do genoma humano para reduzir o tempo de alinhamento
	if [ ! -f $HUMANREFMMI ]; then
		minimap2 -d $HUMANREFMMI $HUMANREFSEQ
	fi
	echo -e "\nExecutando minimap2 & samtools para filtrar as reads do genoma humano..."
	[ ! -d "${READSLEVELDIR}" ] && mkdir -vp ${READSLEVELDIR}
	[ ! -d "${ASSEMBLYDIR}" ] && mkdir -vp ${ASSEMBLYDIR}
	# Loop para analisar todos barcodes, um de cada vez
	for i in $(find ${PRINSEQDIR} -type f -name "*.good.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		echo -e "\nCarregando os dados ${i}..."
	    	# Alinha as reads contra o arquivo indice do genoma humano e ordena os segmentos
	    	minimap2 -ax map-ont -t ${THREADS} ${HUMANREFMMI} ${PRINSEQDIR}/${i}.good.fastq | samtools sort -@ ${THREADS} -o ${READSLEVELDIR}/${i}.sorted.bam -
	    	# Indexa o arquivo para acesso mais rápido
	    	samtools index -@ ${THREADS} ${READSLEVELDIR}/${i}.sorted.bam
	    	# Filtra os segmentos não mapeados Flag 4 (-f 4) para um novo arquivo unmapped.sam 
	    	samtools view -bS -f 4 ${READSLEVELDIR}/${i}.sorted.bam > ${READSLEVELDIR}/${i}.unmapped.bam -@ ${THREADS}
		# Salva os dados no formato .fastq
		samtools fastq ${READSLEVELDIR}/${i}.unmapped.bam > ${READSLEVELDIR}/${i}.unmapped.fastq -@ ${THREADS}
	done
}

function autocorrection () {
	# Autocorreção das reads
	$PRINSEQDIR=$1
	$THREADS=$2
	$READSLEVELDIR=$3
	echo -e "\nExecutando minimap2 & racon para autocorreção das reads contra a sequencia consenso..."
	for i in $(find ${PRINSEQDIR} -type f -name "*.good.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		echo -e "\nCarregando os dados ${i}..."
		# Alinhar todas as reads com elas mesmas para produzir sequencias consenso a partir do overlap de reads
		minimap2 -ax ava-ont -t ${THREADS} ${READSLEVELDIR}/${i}.unmapped.fastq ${READSLEVELDIR}/${i}.unmapped.fastq > ${READSLEVELDIR}/${i}.overlap.sam
		# Correção de erros a partir das sequencias consenso
	 	racon -t ${THREADS} -f -u ${READSLEVELDIR}/${i}.unmapped.fastq ${READSLEVELDIR}/${i}.overlap.sam ${READSLEVELDIR}/${i}.unmapped.fastq > ${READSLEVELDIR}/${i}.corrected.fasta
	done
}

function kraken () {
	# Classificação taxonômica utilizando Kraken2
	$KRAKENDB=$1
	$THREADS=$2
	$READSLEVELDIR=$3
	$RUNNAME=$4
	$MODEL=$5
	echo -e "\nExecutando o Kraken2..."
	for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		# kraken2 --db ${KRAKENDB} --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --report-minimizer-data --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta
		echo -e "\nCarregando os dados ${i}..."
		kraken2 --db ${KRAKENDB} --quick --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta
		echo -e "\nResultados ${i}"
		~/scripts/kraken2_quick_report.sh "${RUNNAME}_${MODEL}" "${i}"
	done
}

# Define as etapas e argumentos de cada workflow
workflowList=(
	'sequencing_summary1:RAWDIR basecalling:RAWDIR;BASECALLDIR;MODEL'
	'sequencing_summary1:RAWDIR basecalling:RAWDIR;BASECALLDIR;MODEL demux_cat1:BASECALLDIR;DEMUXDIR;DEMUXCATDIR sequencing_summary2:RESULTSDIR;RUNNAME;BASECALLDIR;DEMUXDIR;QSCORE;LENGTH primer_removal:CUTADAPTDIR;PRIMER;DEMUXCATDIR qc_filter1:CUTADAPTDIR;LENTGH;NANOFILTERDIR;PRINSEQDIR;QUERYDIR blast:BLASTDIR;QUERYDIR;RUNNAME;MODEL'
	'sequencing_summary1:RAWDIR basecalling:RAWDIR;BASECALLDIR;MODEL demux_cat2:BASECALLDIR;DEMUXDIR;DEMUXCATDIR sequencing_summary2:RESULTSDIR;RUNNAME;BASECALLDIR;DEMUXDIR;QCORE;LENGTH qc_filter2:DEMUXCATDIR;LENGTH;NANOFILTERDIR;PRINSEQDIR human_filter:HUMANREFDIR;READSLEVELDIR;ASSEMBLYDIR;PRINSEQDIR autocorrection:PRINSEQDIR;THREADS;READSLEVELDIR kraken:KRAKENDB;THREADS;READSLEVELDIR;RUNNAME;MODEL'
)

# Índice do array 0..n
indice=$(expr $WF - 1)

# Executa as etapas do workflow selecionado
echo "Executando o workflow WF$WF..."
echo "Passos do WF$WF: ${workflowList[$indice]}"
echo "Libray: $RUNNAME"
echo "Modelo: $MODEL"
# Separa cada etapa do workflow no vetor steps
read -r -a steps <<< "${workflowList[0]}"
for call_func in "${steps[@]}"; do 
	echo "Step: $i"
	# Separo o nome da função dos argumentos
	IFS=":" read -r -a func_name <<< $call_func
	args=$(echo "${func_name[1]}" | tr ";" " ")
	# Obtem os valores de cada argumento
	args_values=""
	for j in $args; do
		echo "Argumentos separados: $j - Valor do argumento: ${!j}"
		args_values="$args_values ${!j}"
	done
	echo "Executando a função ${func_name[0]}"..."
	echo "Argumentos passados: $args...
	echo "Valores dos argumentos: $args_values..."
	# Executa o código e estima o tempo de execução
	export -f "$call_func"
	echo "$call_fun $args_values" | /usr/bin/time -o ~/performance-analysis/${RUNNAME}_${i}.time /bin/bash
done

exit 0
