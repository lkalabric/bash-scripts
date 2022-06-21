#!/bin/bash

# script: metagenomics4.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 09 JUN 2022
# última atualização: 09 JUN 2022
# versão 1: modulariza as etapas do workflow e permite criar diferentes wokflows executado cada etapa como uma função

# Descrição de cada etapa disponível para construção dos workflows
# sequencing_summary1
# sequencing_summary2
# basecalling
# demux
# demux_headcrop
# primer_removal
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
RAWDIR="${HOME}/data/${RUNNAME}"
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
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODEL}"
[ ! -d "${RESULTSDIR}" ] && mkdir -vp ${RESULTSDIR}
BASECALLDIR="${RESULTSDIR}/BASECALL"
# Remove as pastas de resultados anteriores antes das análises
[ ! -d "${RESULTSDIR}/wf${WF}" ] && rm -r "${RESULTSDIR}/wf${WF}"
DEMUXDIR="${RESULTSDIR}/wf${WF}/DEMUX"
DEMUXCATDIR="${RESULTSDIR}/wf${WF}/DEMUX_CAT"
CUTADAPTDIR="${RESULTSDIR}/wf${WF}/CUTADAPT"
NANOFILTDIR="${RESULTSDIR}/wf${WF}/NANOFILT"
PRINSEQDIR="${RESULTSDIR}/wf${WF}/PRINSEQ"
QUERYDIR="${RESULTSDIR}/wf${WF}/QUERY"
BLASTDIR="${RESULTSDIR}/wf${WF}/BLAST"
READSLEVELDIR="${RESULTSDIR}/wf${WF}/READS_LEVEL"
CONTIGLEVELDIR="${RESULTSDIR}/wf${WF}/CONTIGS_LEVEL"
ASSEMBLYDIR="${RESULTSDIR}/wf${WF}/ASSEMBLY"

# Pausa a execução para debug
# read -p "Press [Enter] key to continue..."

# Parâmetros de qualidade mínima
QSCORE=9	# Default Fast min_qscore=8; Hac min_qscore=9; Sup min_qscore=10
LENGTH=100

# Parâmetro de otimização minimap2, samtools, racon e kraken2
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros minimap2 
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -P ${HUMANREFDIR}
HUMANREFSEQ="${HUMANREFDIR}/GRCh38.p13.genome.fa.gz"
HUMANREFMMI="${HUMANREFDIR}/GRCh38.p13.genome.mmi"
# Cria o arquivo índice do genoma humano para reduzir o tempo de alinhamento
if [ ! -f $HUMANREFMMI ]; then
	minimap2 -d $HUMANREFMMI $HUMANREFSEQ
fi

function sequencing_summary1 () {
	# Sumario do sequenciamento (dados disponíveis no arquivo report*.pdf)
	echo "Sumário da corrida"
	echo "Total files:"
	ls $(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;) | wc -l
	echo "Total reads:"
	# h5ls "$(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;)"/*.fast5 | wc -l
}

function basecalling () {
	# Basecalling (comum a todos workflows)
	# Parâmetros Guppy basecaller (ONT)
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
	# Comando para pycoQC version 2.5
	if [ ! -f "${RESULTSDIR}/basecalling_pycoqc.html" ]; then
		echo -e "\nExecutando pycoQC no sequencing summary com o parâmetro default QSCORE=8..."
		pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -o "${RESULTSDIR}/basecalling_wf${wf}_pycoqc.html" --report_title $RUNNAME --min_pass_qual ${QSCORE}
	fi
	if [ ! -f "${RESULTSDIR}/barcoding_pycoqc.html" ]; then
		echo -e "\nExecutando pycoQC no sequencing e barecoder summaries utilizandos os LENGHT=100 e QSCORE=9..."
		pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/barcoding_summary.txt" -o "${RESULTSDIR}/barcoding_wf${wf}_pycoqc.html" --report_title $RUNNAME --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}
	fi
}

function primer_removal () {
	# Remoção dos primers
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
	# Link: https://chipster.csc.fi/manual/prinseq-complexity-filter.html
	[ ! -d ${PRINSEQDIR} ] && mkdir -vp ${PRINSEQDIR}
	echo -e "\nExecutando prinseq-lite.pl..."
	for i in $(find ${NANOFILTDIR} -type f -exec basename {} .fastq \; | sort); do
		echo -e "\nResultados ${i}..."
		prinseq-lite.pl -fastq "${NANOFILTDIR}/${i}.fastq" -out_good "${PRINSEQDIR}/${i}.good" -out_bad "${PRINSEQDIR}/${i}.bad -graph_data" "${PRINSEQDIR}/${i}.gd" -no_qual_header -lc_method dust -lc_threshold 40
		# Resultados disponíveis no report do Prinseq (Good sequences)
	done
}

function blast () {
	# Classificação taxonômica utilizando blastn
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
	[ ! -d "${READSLEVELDIR}" ] && mkdir -vp ${READSLEVELDIR}
	[ ! -d "${ASSEMBLYDIR}" ] && mkdir -vp ${ASSEMBLYDIR}
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
	echo -e "\nExecutando o Kraken2..."
	for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		# kraken2 --db ${KRAKENDB} --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --report-minimizer-data --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta
		echo -e "\nCarregando os dados ${i}..."
		kraken2 --db ${KRAKENDB} --quick --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta
		echo -e "\nResultados ${i}"
		~/scripts/kraken2_quick_report.sh "${RUNNAME}_${MODEL}" "${i}"
	done
}

# Define as etapas de cada workflow
workflowList=(
	'sequencing_summary1 basecalling'
	'sequencing_summary1 basecalling demux sequencing_summary2 primer_removal qc_filter1 qc_filter2 blast'
	'sequencing_summary1 basecalling demux_headcrop sequencing_summary2 qc_filter1 qc_filter2 human_filter autocorrection kraken'
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
# Separa os passos do workflow
read -r -a steps <<< "${workflowList[$indice]}"

for call_func in ${steps[@]}; do
	echo "Executando a função $call_func..."
	eval $call_func	
done

exit 0

# Análise de cobertura do sequenciamento
# Genomas referência e plot de cobertura
CHIKVREFSEQ="${REFSEQ}/Togaviridae/NC_004162.2_CHIKV-S27.fasta"
DENV1REFSEQ="${REFSEQ}/Flaviviridae/NC_001477.1_DENV1.fasta"
DENV2REFSEQ="${REFSEQ}/Flaviviridae/NC_001474.2_DENV2.fasta"
DENV3REFSEQ="${REFSEQ}/Flaviviridae/NC_001475.2_DENV3.fasta"
DENV4REFSEQ="${REFSEQ}/Flaviviridae/NC_002640.1_DENV4.fasta"
ZIKVREFSEQ="${REFSEQ}/Flaviviridae/NC_012532.1_ZIKV.fasta"
HIV1REFSEQ="${REFSEQ}/Retroviridae/NC_001802.1_HIV1.fasta"

# Mapeamento CHIKV
#for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${CHIKVREFSEQ} ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.chikv.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.chikv.sorted.bam > ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${CHIKVREFSEQ} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.chikv -n N -i ${i}
#done

# Mapeamento DENV1
#for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV1REFSEQ} ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv1.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv1.sorted.bam > ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV1REFSEQ} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv1 -n N -i ${i}
#done

# Mapeamento DENV2
#for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV2REFSEQ} ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv2.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv2.sorted.bam > ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV2REFSEQ} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv2 -n N -i ${i}
#done

# Mapeamento DENV3
#for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV3REFSEQ} ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv3.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv3.sorted.bam > ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV3REFSEQ} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv3 -n N -i ${i}
#done

# Mapeamento DENV4
#for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV4REFSEQ} ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv4.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv4.sorted.bam > ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV4REFSEQ} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv4 -n N -i ${i}
#done

# Mapeamento ZIKV
#for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${ZIKVREFSEQ} ${READSLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.zikv.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.zikv.sorted.bam > ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${ZIKVREFSEQ} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.zikv -n N -i ${i}
#done

#source activate plot
#for i in $(find ${ASSEMBLYDIR} -type f -name "*.sorted.mapped.bam" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	fastcov ${ASSEMBLYDIR}/barcode*.chikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_chikv.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.chikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_chikv_log.pdf
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
#done

#
# Contig-level taxid
#


