#!/bin/bash

# script: metagenomics5.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 09 JUN 2022
# última atualização: 27 JUN 2022
# versão 5: modulariza as etapas do workflow e permite criar diferentes wokflows executado cada etapa como uma função e analisar o tempo de execução de cada etapa

# Descrição de cada etapa disponível para construção dos workflows
# sequencing_summary1
# sequencing_summary2
# basecalling
# demux_cat1
# demux_cat2
# primer_removal
# qc_filter1
# qc_filter2
# human_filter
# autocorrection
# blastn_local
# kraken_local

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=$3		# Workflow de bioinformatica 1, 2 ou 3
if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, número do worflow ou modelo Guppy Basecaller!"
	echo "Sintáxe: ./metagenomics4time.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,3>"
	exit 0
fi

# Caminho de INPUT dos dados fast5
RAWDIR="${HOME}/data/${RUNNAME}"
if [ ! -d $RAWDIR ]; then
	echo "Pasta de dados não encontrada!"
	exit 0
fi

# Caminho de INPUT dos bancos de dados
HUMANREFDIR=${HOME}/data/GRCh38
REFSEQDIR=${HOME}/data/REFSEQ
BLASTDBDIR="${HOME}/data/BLAST_DB"
KRAKENDBDIR="${HOME}/data/KRAKEN2_DB" # Substituir pelo nosso banco de dados se necessário KRAKEN2_USER_DB

# Parâmetros minimap2 
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -P ${HUMANREFDIR}
HUMANREFSEQ="${HUMANREFDIR}/GRCh38.p13.genome.fa.gz"
HUMANREFMMI="${HUMANREFDIR}/GRCh38.p13.genome.mmi"
# Cria o arquivo índice do genoma humano para reduzir o tempo de alinhamento
if [ ! -f $HUMANREFMMI ]; then
	minimap2 -d $HUMANREFMMI $HUMANREFSEQ
fi

# Caminhos de OUTPUT das análises
echo "Preparando pastas para (re-)análise dos dados..."
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}/time"
[ ! -d "${RESULTSDIR}" ] && mkdir -vp "${RESULTSDIR}"
# Reseta a pasta de resultados anteriores da worflow 
[ -d "${RESULTSDIR}/wf${WF}" ] && rm -r "${RESULTSDIR}/wf${WF}"; mkdir -vp "${RESULTSDIR}/wf${WF}"

# Pausa a execução para debug
# read -p "Press [Enter] key to continue..."

function sequencing_summary1 () {
  RAWDIR=$1
  # Sumario do sequenciamento (dados disponíveis no arquivo report*.pdf)
  echo "Sumário da corrida"
  echo "Total files:"
  ls $(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;) | wc -l
  echo "Total reads:"
  # h5ls "$(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;)"/*.fast5 | wc -l
}

function basecalling () {
	# Basecalling  (comum a todos workflows)
	RAWDIR=$1
	MODEL=$2
	RESULTSDIR=$3
	BASECALLDIR="${RESULTSDIR}/BASECALL"
	# Parâmetros Guppy basecaller (ONT)
	QSCORE=9 	# Defalut Fast min_qscore=8; Hac min_qscore=9; Sup min_qscore=10
	CONFIG="dna_r9.4.1_450bps_${MODEL}.cfg" #dna_r9.4.1_450bps_fast.cfg dna_r9.4.1_450bps_hac.cfg dna_r9.4.1_450bps_sup.cfg
	# Parâmetros para otimização do Guppy basecaller (ONT) obtidos pelo benckmark utilizadando o LAPTOP-Yale
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
	if [ ! -d $BASECALLDIR ]; then
		mkdir -vp $BASECALLDIR
		echo -e "\nExecutando guppy_basecaller..."
		# Comando para guppy_basecaller usando GPU
		guppy_basecaller -r -i ${RAWDIR} -s "${BASECALLDIR}" -c ${CONFIG} -x auto --min_qscore ${QSCORE} --gpu_runners_per_device ${GPUPERDEVICE} --chunk_size ${CHUNCKSIZE} --chunks_per_runner ${CHUNKPERRUNNER} --verbose_logs
	fi
}

# Para debug
# if false; then # Desvio para execução rápida
# fi # Fim do desvio para execução rápida

function demux () {
	# Demultiplex, adapter removal & sem headcrop 18 para uso do cutadapt
	RESULTSDIR=$1
	WF=$2
	BASECALLDIR="${RESULTSDIR}/BASECALL"
	DEMUXDIR="${RESULTSDIR}/wf${WF}/DEMUX"
	DEMUXCATDIR="${RESULTSDIR}/wf${WF}/DEMUX_CAT"
	# Parâmetros Guppy barcoder (ONT)
	ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
	if [ ! -d $DEMUXDIR ]; then
		mkdir -vp $DEMUXDIR
		echo -e "\nExecutando guppy_barcoder..."
		guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --trim_barcodes  
		# Renomeia a pasta contendo as reads unclassified para barcode00 para análise
		[ -d "${DEMUXDIR}/unclassified" ] && mv "${DEMUXDIR}/unclassified" "${DEMUXDIR}/barcode00"
		# Concatena todos arquivos .fastq de cada barcode em um arquivo .fastq único
		[ ! -d ${DEMUXCATDIR} ] && mkdir -vp ${DEMUXCATDIR}
		for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
			[ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > "${DEMUXCATDIR}/${i}.fastq"
		done
	fi
}

function demux_headcrop () {
	# Demultiplex, adapter removal com headcrop 18 sem uso do cutadapt
	RESULTSDIR=$1
	WF=$2
	BASECALLDIR="${RESULTSDIR}/BASECALL"
	DEMUXDIR="${RESULTSDIR}/wf${WF}/DEMUX"
	DEMUXCATDIR="${RESULTSDIR}/wf${WF}/DEMUX_CAT"
	# Parâmetros Guppy barcoder (ONT)
	TRIMADAPTER=18
	ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
	if [ ! -d $DEMUXDIR ]; then
		mkdir -vp $DEMUXDIR
		echo -e "\nExecutando guppy_barcoder..."
		guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER}
		# Renomeia a pasta contendo as reads unclassified para barcode00 para análise
		[ -d "${DEMUXDIR}/unclassified" ] && mv "${DEMUXDIR}/unclassified" "${DEMUXDIR}/barcode00"
		# Concatena todos arquivos .fastq de cada barcode em um arquivo .fastq único
		[ ! -d ${DEMUXCATDIR} ] && mkdir -vp ${DEMUXCATDIR}
		for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
			[ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > "${DEMUXCATDIR}/${i}.fastq"
		done
	fi
}

function sequencing_summary2 () {
	# pycoQC summary
	$RESULTSDIR=$1
	WF=$2
	BASECALLDIR="${RESULTSDIR}/BASECALL"
	DEMUXDIR="${RESULTSDIR}/wf${WF}/DEMUX"
	# Parâmetros de qualidade mínima
	QSCORE=9	# Defalut Fast min_qscore=8; Hac min_qscore=9; Sup min_qscore=10
	LENGTH=100
	# Comando para pycoQC version 2.5
	if [ ! -f "${RESULTSDIR}/basecalling_wf${WF}_pycoqc.html" ]; then
		echo -e "\nExecutando pycoQC no sequencing summary com o parâmetro default QSCORE=8..."
		pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -o "${RESULTSDIR}/basecalling_wf${WF}_pycoqc.html" --report_title ${RESULTSDIR} --min_pass_qual ${QSCORE}
	fi
	if [ ! -f "${RESULTSDIR}/barcoding_wf${WF}_pycoqc.html" ]; then
		echo -e "\nExecutando pycoQC no sequencing e barecoder summaries utilizandos os LENGHT=100 e QSCORE=9..."
		pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/barcoding_summary.txt" -o "${RESULTSDIR}/barcoding_wf${WF}_pycoqc.html" --report_title ${RESULTSDIR} --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}
	fi
}

function primer_removal () {
	# Remoção dos primers
	RESULTSDIR=$1
	WF=$2
	DEMUXCATDIR="${RESULTSDIR}/wf${WF}/DEMUX_CAT"
	CUTADAPTDIR="${RESULTSDIR}/wf${WF}/CUTADAPT"
	PRIMER="GTTTCCCACTGGAGGATA"
	[ ! -d ${CUTADAPTDIR} ] && mkdir -vp ${CUTADAPTDIR}
	echo -e "\nExecutando cutadapt..."
	for i in $(find ${DEMUXCATDIR} -type f -exec basename {} .fastq \; | sort); do
		cutadapt -g ${PRIMER} -e 0.2 --discard-untrimmed -o "${CUTADAPTDIR}/${i}.fastq" "${DEMUXCATDIR}/${i}.fastq"
		echo -e "\nResultados ${i} $(grep -c "runid" ${CUTADAPTDIR}/${i}.fastq | cut -d : -f 2 | awk '{s+=$1} END {printf "%.0f\n",s}')"
	done
}

function qc_filter1 () {
	# Filtro por tamanho
	RESULTSDIR=$1
	WF=$2
	CUTADAPTDIR="${RESULTSDIR}/wf${WF}/CUTADAPT"
	NANOFILTDIR="${RESULTSDIR}/wf${WF}/NANOFILT"
	# Parâmetros de qualidade mínima
	LENGTH=100
	source activate ngs
	[ ! -d ${NANOFILTDIR} ] && mkdir -vp ${NANOFILTDIR}
	echo -e "\nExecutando NanoFilt..."
	case $WF in
	  2)
		for i in $(find "${CUTADAPTDIR}" -type f -exec basename {} .fastq \; | sort); do
			NanoFilt -l ${LENGTH} < "${CUTADAPTDIR}/${i}.fastq" > "${NANOFILTDIR}/${i}.fastq" 
		done
	  ;;
	  3)
		for i in $(find "${DEMUXCATDIR}" -type f -exec basename {} .fastq \; | sort); do
			NanoFilt -l ${LENGTH} < "${DEMUXCATDIR}/${i}.fastq" > "${NANOFILTDIR}/${i}.fastq" 
		done	  ;;
	  *)
		echo "Erro na função qc_filter1!"
		exit 0;
	  ;;
	esac
}

function qc_filter2 () {
	# Filtro de complexidade
	RESULTSDIR=$1
	WF=$2
	NANOFILTDIR="${RESULTSDIR}/wf${WF}/NANOFILT"
	PRINSEQDIR="${RESULTSDIR}/wf${WF}/PRINSEQDIR"
	# Link: https://chipster.csc.fi/manual/prinseq-complexity-filter.html
	[ ! -d ${PRINSEQDIR} ] && mkdir -vp ${PRINSEQDIR}
	echo -e "\nExecutando prinseq-lite.pl..."
	for i in $(find ${NANOFILTDIR} -type f -exec basename {} .fastq \; | sort); do
		echo -e "\nResultados ${i}..."
		prinseq-lite.pl -fastq "${NANOFILTDIR}/${i}.fastq" -out_good "${PRINSEQDIR}/${i}.good" -out_bad "${PRINSEQDIR}/${i}.bad -graph_data" "${PRINSEQDIR}/${i}.gd" -no_qual_header -lc_method dust -lc_threshold 40
		# Resultados disponíveis no report do Prinseq (Good sequences)
	done
}

function blastn_local () {
	# Classificação taxonômica utilizando blastn
	RESULTSDIR=$1
	WF=$2
	PRINSEQDIR="${RESULTSDIR}/wf${WF}/PRINSEQDIR"
	QUERYDIR="${RESULTSDIR}/wf${WF}/QUERY"
	BLASTDIR="${RESULTSDIR}/wf${WF}/BLAST"
	# Preparação do BLASTDB local
	# Script: makeblastdb_refseq.sh
		# Concatena todas as REFSEQs num arquivo refseq.fasta único e cria o BLASTDB
		# Extrai do arquvio refseq.fasta a lista acesso refseq.acc
		# Cria a partir do arquivo refseq.acc o arquivo refseq.map que mapeia os taxid (números que identificam as espécies taxonômica)
	# Converte arquivos .fastq em .fasta para query no blastn
	[ ! -d "${QUERYDIR}" ] && mkdir -vp ${QUERYDIR}
	for i in $(find "${PRINSEQDIR}"/*.good.fastq -type f -exec basename {} .good.fastq \; | sort); do
		sed -n '1~4s/^@/>/p;2~4p' "${PRINSEQDIR}/${i}.good.fastq" > "${QUERYDIR}/${i}.fasta"
	done
	# Busca as QUERIES no BLASTDB local
	[ ! -d ${BLASTDIR} ] && mkdir -vp ${BLASTDIR}
	echo -e "\nExecutando blastn..."
	for i in $(find ${QUERYDIR} -type f -exec basename {} .fasta \; | sort); do
		echo -e "\nAnalisando dados ${BLASTDIR}/${i}..."
		blastn -db "${BLASTDBDIR}/refseq" -query "${QUERYDIR}/${i}.fasta" -out "${BLASTDIR}/${i}.blastn" -outfmt "6 sacc staxid" -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
		# Busca remota
		# blastn -db nt -remote -query ${QUERYDIR}/${i}.fasta -out ${BLASTDIR}/${i}.blastn -outfmt "6 qacc saccver pident sscinames length mismatch gapopen evalue bitscore"  -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
		echo -e "\nResultados ${i}"
		~/scripts/blast_report.sh "${BLASTDIR}/${i}.blastn"
	done
}

function human_filter () {
	# Remoção das reads do genoma humano
	RESULTSDIR=$1
	WF=$2
	$HUMANREFDIR=$3
	PRINSEQDIR="${RESULTSDIR}/wf${WF}/PRINSEQDIR"
	READSLEVELDIR="${RESULTSDIR}/wf${WF}/READS_LEVEL"
	[ ! -d "${READSLEVELDIR}" ] && mkdir -vp ${READSLEVELDIR}
	# Parâmetro de otimização minimap2, samtools, racon e kraken2
	THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"
	echo -e "\nExecutando minimap2 & samtools para filtrar as reads do genoma humano..."
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
	RESULTSDIR=$1
	WF=$2
	PRINSEQDIR="${RESULTSDIR}/wf${WF}/PRINSEQDIR"
	READSLEVELDIR="${RESULTSDIR}/wf${WF}/READS_LEVEL"
	# Parâmetro de otimização minimap2, samtools, racon e kraken2
	THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"
	echo -e "\nExecutando minimap2 & racon para autocorreção das reads contra a sequencia consenso..."
	for i in $(find ${PRINSEQDIR} -type f -name "*.good.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		echo -e "\nCarregando os dados ${i}..."
		# Alinhar todas as reads com elas mesmas para produzir sequencias consenso a partir do overlap de reads
		minimap2 -ax ava-ont -t ${THREADS} ${READSLEVELDIR}/${i}.unmapped.fastq ${READSLEVELDIR}/${i}.unmapped.fastq > ${READSLEVELDIR}/${i}.overlap.sam
		# Correção de erros a partir das sequencias consenso
	 	racon -t ${THREADS} -f -u ${READSLEVELDIR}/${i}.unmapped.fastq ${READSLEVELDIR}/${i}.overlap.sam ${READSLEVELDIR}/${i}.unmapped.fastq > ${READSLEVELDIR}/${i}.corrected.fasta
	done
}

function kraken_local () {
	# Classificação taxonômica utilizando Kraken2
	RESULTSDIR=$1
	WF=$2
	KRAKENDBDIR=$3
	READSLEVELDIR="${RESULTSDIR}/wf${WF}/READS_LEVEL"
	# Parâmetro de otimização minimap2, samtools, racon e kraken2
	THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"
	echo -e "\nExecutando o Kraken2..."
	for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		# kraken2 --db ${KRAKENDBDIR} --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --report-minimizer-data --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta
		echo -e "\nCarregando os dados ${i}..."
		kraken2 --db ${KRAKENDBDIR} --quick --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta
		echo -e "\nResultados ${i}"
		~/scripts/kraken2_quick_report.sh "${READSLEVELDIR}/${i}_report.txt"
	done
}

function assembly () {
# Faz a análise de cobertura e montagem das reads em sequencias referências
	CONTIGLEVELDIR="${RESULTSDIR}/wf${WF}/CONTIGS_LEVEL"
	ASSEMBLYDIR="${RESULTSDIR}/wf${WF}/ASSEMBLY"
	[ ! -d "${ASSEMBLYDIR}" ] && mkdir -vp ${ASSEMBLYDIR}
	case $WF in
		2)
	  		echo "Montando as sequencias do WF $WF..."		
	  	;;
	  	3)
	    		echo "Montando as sequencias do WF $WF..."
	    	 ;;
	  	*)
	    		echo "Dados não disponíveis para o WF $WF!"
			exit 1
		;;
	esac
}

# Define as etapas e argumentos de cada workflow
workflowList=(
	'sequencing_summary1:RAWDIR basecalling:RAWDIR;MODEL;RESULTSDIR'
	'sequencing_summary1:RAWDIR basecalling:RAWDIR;MODEL;RESULTSDIR demux:RESULTSDIR;WF sequencing_summary2:RESULTSDIR;WF primer_removal:RESULTSDIR;WF qc_filter1:RESULTSDIR;WF qc_filter2:RESULTSDIR;WF blastn_local:RESULTSDIR;WF'
	'sequencing_summary1:RAWDIR basecalling:RAWDIR;MODEL;RESULTSDIR demux_headcrop:RESULTSDIR;WF sequencing_summary2:RESULTSDIR;WF qc_filter1:RESULTSDIR;WF qc_filter2:RESULTSDIR;WF human_filter:RESULTSDIR;WF;HUMANREFDIR autocorrection:RESULTSDIR;WF kraken_local:RESULTSDIR;WF;KRAKENDBDIR'
)


# Executa as etapas do workflow selecionado

# Validação do WF
if [[ $WF -gt ${#workflowList[@]} ]]; then
	echo "Workflow não definido!"
	exit 0;
fi
# Índice para o array workflowList 0..n
indice=$(expr $WF - 1)

# Execução das análises propriamente ditas a partir do workflow selecionado
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
