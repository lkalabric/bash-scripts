#!/bin/bash

# script: metagenomics_wf3_v1_lk.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 25 AGO 2021
# versão 3.1: corrigiu o wf 3

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODELO=$2	# Modelo de basecalling fast hac sup
WF=$3		# Workflow de bioinformatica 1, 2 ou 31

if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados, número do worflow ou modelo Guppy Basecaller!"
	echo "Sintáxe: ./metagenomic_v3_lk.sh <LIBRARY> <MODELO:fast,hac,sup> <WF: 1,2,31>"
	exit 0
fi
# Caminho de INPUT dos dados fast5
RAWDIR="${HOME}/data/${RUNNAME}" # Se análise começar com o Barcoder
 
if [ ! -d $RAWDIR ]; then
	echo "Pasta de dados não encontrada!"
	exit 0
fi

# Parâmetros Guppy (ONT)
CONFIG="dna_r9.4.1_450bps_${MODELO}.cfg" #dna_r9.4.1_450bps_fast.cfg dna_r9.4.1_450bps_sup.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"

# Caminhos de output das análises
RESULTSDIR="${HOME}/ngs-analysis/${RUNNAME}_${MODELO}"
rm -r ${RESULTSDIR}
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

# Caminho para os bancos de dados
REFSEQDIR=${HOME}/data/REFSEQ
HUMANREFDIR=${HOME}/data/GRCh38
BLASTDBDIR="${HOME}/data/BLAST_DB"
KRAKENDB="${HOME}/data/KRAKEN2_DB" # Substituir pelo nosso banco de dados se necessário KRAKEN2_USER_DB
 
# Parâmetros para otimização do Guppy basecaller para modelo fast (benckmark)
GPUPERDEVICE=4		# Parâmetros otimizados a partir do benchmark GPU LAPTOP-Yale
CHUNCKSIZE=1000		# Parâmetros otimizados a partir do benchmark GPU LAPTOP-Yale
CHUNKPERRUNNER=50	# Parâmetros otimizados a partir do benchmark GPU LAPTOP-Yale

# Parâmetro de otimização minimap2, samtools, racon e kraken2
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros Barcoder ou Cutadapt para remoção do primer
TRIMADAPTER=18
PRIMER="GTTTCCCACTGGAGGATA"

# Parâmetros de qualidade mínima
QSCORE=9
LENGTH=100

# Parâmetros minimap2 
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -P ${HUMANREFDIR}
HUMANREFSEQ="${HOME}/data/GRCh38/GRCh38.p13.genome.fa.gz"
HUMANREFMMI="${HOME}/data/GRCh38/GRCh38.p13.genome.mmi"
if [ ! -f $HUMANREFMMI ]; then
	minimap2 -d $HUMANREFMMI $HUMANREFSEQ
fi

# Todos os WFs
# Sumarios dos dados
echo "Sumário dos dados"
ls $(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;) | wc -l
h5ls $(find ${RAWDIR} -type f -name "*.fast5" -exec diname {} \;) | wc -l

# Basecalling
if [ ! -d $BASECALLDIR ]; then
	mkdir -vp $BASECALLDIR
	# Step 1 - Basecalling (comum a todos workflows)
	# Esta etapa está sendo realizada pelo script guppy_gpu_v1_ag.sh no LAPTOP-Yale
	echo Executando guppy_basecaller...
	# Comando para guppy_basecaller usando GPU
	guppy_basecaller -r -i ${RAWDIR} -s "${BASECALLDIR}" -c ${CONFIG} -x auto  --gpu_runners_per_device ${GPUPERDEVICE} --chunk_size ${CHUNCKSIZE} --chunks_per_runner ${CHUNKPERRUNNER} --verbose_logs
	echo "Basecalling concluido com sucesso!"
fi

# WF 1 - Classificação Taxonômica pelo Epi2ME
# Copiar a pasta /pass para o Epi2ME
if [[ $WF -eq 1 ]]; then
	exit 1
fi

# if false; then # Descio para execução rápida

# Workflow 2 & 31
# Step 2 - Demultiplex & adapter removal (comum a todos workflows)
if [ ! -d $DEMUXDIR ]; then
	mkdir -vp $DEMUXDIR
	echo "Executando guppy_barcoder..."
	#bm1 guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --trim_barcodes
	#bm2 guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --detect_mid_strand_adapter --trim_barcodes
	#bm3 guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --trim_barcodes
	#bm4 guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --detect_mid_strand_barcodes --trim_barcodes  
	#bm6 guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --barcode_kits EXP-NBD104 --require_barcodes_both_ends --detect_mid_strand_barcodes --trim_barcodes  
	if [[ $WF -eq 2 ]];
	then
		#bm5 WF2 sem headcrop 18 (uso cutadapt)
		guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --detect_mid_strand_adapter --trim_barcodes  
	else
		#bm7 WF31 com headcrop 18
		guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER}
	fi
fi
# Move a pasta contendo as reads unclassified para barcode00
[ -d "${DEMUXDIR}/unclassified" ] && mv "${DEMUXDIR}/unclassified" "${DEMUXDIR}/barcode00"

# fi # Fim do desvio para execução rápida

# Step 3 - QC
#echo "Executando pycoQC..."
source activate ngs
if [ ! -f "${RESULTSDIR}/${RUNNAME}_pycoqc.html" ]; then
	# Comando para pycoQC version 2.5
	pycoQC -q -f "${BASECALLDIR}/sequencing_summary.txt" -b "${DEMUXDIR}/barcoding_summary.txt" -o "${RESULTSDIR}/${RUNNAME}_pycoqc.html" --report_title $RUNNAME --min_pass_qual ${QSCORE} --min_pass_len ${LENGTH}
fi

# Concatena todos arquivos .fastq de cada barcode em um arquivo .fastq único
[ ! -d ${DEMUXCATDIR} ] && mkdir -vp ${DEMUXCATDIR}
for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
    [ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > "${DEMUXCATDIR}/${i}.fastq"
done

# Workflow 2 apenas
if [[ $WF -eq 2 ]]; then 
	# Step 4 - Remoção dos primers
	echo "Executando cutadapt..."
	[ ! -d ${CUTADAPTDIR} ] && mkdir -vp ${CUTADAPTDIR}
	for i in $(find ${DEMUXCATDIR} -type f -exec basename {} .fastq \;); do
		cutadapt -g ${PRIMER} -e 0.2 --discard-untrimmed -o "${CUTADAPTDIR}/${i}.fastq" "${DEMUXCATDIR}/${i}.fastq"
	done;

	# Step 5 - Filtro por tamanho
	echo "Executando NanoFilt..."
	[ ! -d ${NANOFILTDIR} ] && mkdir -vp ${NANOFILTDIR}
	for i in $(find "${CUTADAPTDIR}" -type f -exec basename {} .fastq \;); do
		NanoFilt -l ${LENGTH} < "${CUTADAPTDIR}/${i}.fastq" > "${NANOFILTDIR}/${i}.fastq" 
	done;

	# Step 6 - Filtro de complexidade
	# Link: https://chipster.csc.fi/manual/prinseq-complexity-filter.html
	echo "Executando prinseq-lite.pl..."
	[ ! -d ${PRINSEQDIR} ] && mkdir -vp ${PRINSEQDIR}
	for i in $(find ${NANOFILTDIR} -type f -exec basename {} .fastq \;); do
		prinseq-lite.pl -fastq "${NANOFILTDIR}/${i}.fastq" -out_good "${PRINSEQDIR}/${i}.good" -out_bad "${PRINSEQDIR}/${i}.bad -graph_data" "${PRINSEQDIR}/${i}.gd" -no_qual_header -lc_method dust -lc_threshold 40
	done;

	# Converte arquivos .fastq em .fasta para query no blastn
	[ ! -d "${QUERYDIR}" ] && mkdir -vp ${QUERYDIR}
	for i in $(find "${PRINSEQDIR}"/*.good.fastq -type f -exec basename {} .good.fastq \;); do
		sed -n '1~4s/^@/>/p;2~4p' "${PRINSEQDIR}/${i}.good.fastq" > "${QUERYDIR}/${i}.fasta"
	done;
		
	# Step 7 - Classificação taxonômica
	# Script: makeblastdb_refseq.sh
	# Concatena todas as REFSEQs num arquivo refseq.fasta único e cria o BLASTDB

	# Script: fasta2acc.sh 
	# Extrai do arquvio refseq.fasta a lista acesso refseq_acc.txt

	# Script: acc2taxid.sh 
	# Cria a partir do arquivo refser_acc.txt o arquivo refseq_map.txt que mapeia os taxid (números que identificam as espécies taxonômica)

	# Busca as reads no BLASTDB local
	echo "Executando blastn..."
	[ ! -d ${BLASTDIR} ] && mkdir -vp ${BLASTDIR}
	for i in $(find ${QUERYDIR} -type f -exec basename {} .fasta \;); do
			echo "Carregando os dados ${BLASTDIR}/${i}..."
			blastn -db "${BLASTDBDIR}/refseq" -query "${QUERYDIR}/${i}.fasta" -out "${BLASTDIR}/${i}.blastn" -outfmt "6 sacc staxid" -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
			# Busca remota
			# blastn -db nt -remote -query ${QUERYDIR}/${i}.fasta -out ${BLASTDIR}/${i}.blastn -outfmt "6 qacc saccver pident sscinames length mismatch gapopen evalue bitscore"  -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
	done;
	exit 2
fi

# Workflow 3.1 apenas
if [[ $WF -eq 31 ]]; then 
	source activate ngs
	# Step 4 - Filtro por tamanho
	echo "Executando NanoFilt..."
	[ ! -d ${NANOFILTDIR} ] && mkdir -vp ${NANOFILTDIR}
	for i in $(find ${DEMUXCATDIR} -type f -exec basename {} .fastq \;); do
		NanoFilt -l ${LENGTH} < "${DEMUXCATDIR}/${i}.fastq" > "${NANOFILTDIR}/${i}.fastq" 
	done;

	# Step 5 - Filtro de complexidade
	# Link: https://chipster.csc.fi/manual/prinseq-complexity-filter.html
	echo "Executando prinseq-lite.pl..."
	[ ! -d ${PRINSEQDIR} ] && mkdir -vp ${PRINSEQDIR}
	for i in $(find ${NANOFILTDIR} -type f -exec basename {} .fastq \;); do
		prinseq-lite.pl -fastq "${NANOFILTDIR}/${i}.fastq" -out_good "${PRINSEQDIR}/${i}.good" -out_bad "${PRINSEQDIR}/${i}.bad -graph_data" "${PRINSEQDIR}/${i}.gd" -no_qual_header -lc_method dust -lc_threshold 40
	done;

	# Step 6 - Remoção das reads do genoma humano
	echo "Executando minimap2, samtools & racon..."
	[ ! -d "${READSLEVELDIR}" ] && mkdir -vp ${READSLEVELDIR}
	[ ! -d "${ASSEMBLYDIR}" ] && mkdir -vp ${ASSEMBLYDIR}
	# Cria o arquivo índice do genoma humano para reduzir o tempo de alinhamento
	minimap2 -d ${HUMANREFMMI} ${HUMANREFSEQ}

	# Loop para analisar todos barcodes, um de cada vez
	for i in $(find ${PRINSEQDIR} -type f -name "*.good.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	    	# Alinha as reads contra o arquivo indice do genoma humano e ordena os segmentos
	    	minimap2 -ax map-ont -t ${THREADS} ${HUMANREFMMI} ${PRINSEQDIR}/${i}.good.fastq | samtools sort -@ 12 -o ${READSLEVELDIR}/${i}.sorted.bam -
	    	# Indexa o arquivo para acesso mais rápido
	    	samtools index -@ ${THREADS} ${READSLEVELDIR}/${i}.sorted.bam
	    	# Filtra os segmento não mapeados Flag 4 (-f 4) para um novo arquivo .bam 
	    	samtools view -bS -f 4 ${READSLEVELDIR}/${i}.sorted.bam > ${READSLEVELDIR}/${i}.unmapped.sam -@ 12
		# Salva os dados no formato .fastq
		samtools fastq -f 4 ${READSLEVELDIR}/${i}.unmapped.sam > ${READSLEVELDIR}/${i}.unmapped.fastq -@ 12
		# Alinhar todas as reads com elas mesmas para produzir sequencias consenso a partir do overlap de reads
		minimap2 -ax ava-ont -t ${THREADS} ${READSLEVELDIR}/${i}.unmapped.fastq ${READSLEVELDIR}/${i}.unmapped.fastq > ${READSLEVELDIR}/${i}.overlap.sam
		# Correção de erros a partir das sequencias consenso
	 	racon -t ${THREADS} -f -u ${READSLEVELDIR}/${i}.unmapped.fastq ${READSLEVELDIR}/${i}.overlap.sam ${READSLEVELDIR}/${i}.unmapped.fastq > ${READSLEVELDIR}/${i}.corrected.fasta
	done

	# Step 7 - Classificação taxonômica
	echo "Executando o Kraken2..."
	for i in $(find ${READSLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
		# kraken2 --db ${KRAKENDB} --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --report-minimizer-data --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta
		echo "Carregando os dados ${READSLEVELDIR}/${i}..."
		kraken2 --db ${KRAKENDB} --quick --threads ${THREADS} --report ${READSLEVELDIR}/${i}_report.txt --output ${READSLEVELDIR}/${i}_output.txt ${READSLEVELDIR}/${i}.corrected.fasta

	done
	exit 31
fi

echo "Workflow $WF concluido com sucesso!"
exit 0

# Pós-análise
# Mapeamento contra genomas referência e plot de cobertura
CHIKVREFSEQ="${REFSEQ}/Togaviridae/NC_004162.2_CHIKV-S27.fasta"
DENV1REFSEQ="${REFSEQ}/Flaviviridae/NC_001477.1_DENV1"
DENV2REFSEQ="${REFSEQ}/Flaviviridae/NC_001474.2_DENV2"
DENV3REFSEQ="${REFSEQ}/Flaviviridae/NC_001475.2_DENV3"
DENV4REFSEQ="${REFSEQ}/Flaviviridae/NC_002640.1_DENV4"
ZIKVREFSEQ="${REFSEQ}/Flaviviridae/NC_012532.1_ZIKV"

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
