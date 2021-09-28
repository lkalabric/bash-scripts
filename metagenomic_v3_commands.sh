#!/bin/bash

# script: metagenomics_wf3_v1_lk.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# data: 25 AGO 2021
# uso: execute cada comando independemente do WF (exceto headcrop 18)

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
if [[ $# -eq 0 ]]; then
	echo "Falta o nome dos dados!"
	exit 0
fi
# Caminho de INPUT dos dados fast5
RAWDIR="${HOME}/data/${RUNNAME}"
 
if [ ! -d $RAWDIR ]; then
	echo "Pasta de dados não encontrada!"
	exit 0
fi

# Workflow de bioinformática (WF=0, apenas comandos)
WF=0

# Ativa os programas do repositório bioconda instalados no ambiente ngs
source activate ngs

# Caminhos de output das análises
BASECALLDIR="${HOME}/ngs-analysis/${RUNNAME}"
mkdir -p "${BASECALLDIR}"
mkdir -p "${BASECALLDIR}/wf${WF}"
DEMUXDIR="${BASECALLDIR}/wf${WF}/DEMUX"
DEMUXCATDIR="${BASECALLDIR}/wf${WF}/DEMUX_CAT"
CUTADAPTDIR="${BASECALLDIR}/wf${WF}/CUTADAPT"
NANOFILTDIR="${BASECALLDIR}/wf${WF}/NANOFILT"
PRINSEQDIR="${BASECALLDIR}/wf${WF}/PRINSEQ"
QUERYDIR="${BASECALLDIR}/wf${WF}/QUERY"
BLASTDIR="${BASECALLDIR}/wf${WF}/BLAST"
READLEVELDIR="${BASECALLDIR}/wf${WF}/LEVEL_READS"
CONTIGLEVELDIR="${BASECALLDIR}/wf${WF}/LEVEL_CONTIGS"
ASSEMBLYDIR="${BASECALLDIR}/wf${WF}/ASSEMBLY"

# Caminho para os bancos de dados
REFSEQDIR=${HOME}/data/REFSEQ
HUMANREFDIR=${HOME}/data/GRCh38
BLASTDBDIR="${HOME}/data/BLAST_DB"
KRAKENDB="${HOME}/data/KRAKEN2_DB" # Substituir pelo nosso banco de dados se necessário KRAKEN2_USER_DB
 
# Parâmetros Guppy (ONT) 
CONFIG="dna_r9.4.1_450bps_hac.cfg" #dna_r9.4.1_450bps_fast.cfg dna_r9.4.1_450bps_sup.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
GPUPERDEVICE=4		# Parâmetros otimizados a partir do benchmark GPU LAPTOP-Yale
CHUNCKSIZE=1000		# Parâmetros otimizados a partir do benchmark GPU LAPTOP-Yale
CHUNKPERRUNNER=50	# Parâmetros otimizados a partir do benchmark GPU LAPTOP-Yale
TRIMADAPTER=18

# Parâmetro de otimização minimap2, samtools, racon e kraken2
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

# Parâmetros de qualidade
QSCORE=9
LENGTH=100

# Parâmetros Cutadapt
PRIMER="GTTTCCCACTGGAGGATA"

# Parâmetros minimap2 
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -P ${HUMANREFDIR}
HUMANREFSEQ="${HOME}/data/GRCh38/GRCh38.p13.genome.fa.gz"
HUMANREFMMI="${HOME}/data/GRCh38/GRCh38.p13.genome.mmi"


if false; then # Desvio para execução a partir daqui

# Comando 1 - guppy_basecaller
# Esta etapa está sendo realizada pelo script guppy_gpu_v1_ag.sh no LAPTOP-Yale
echo Executando guppy_basecaller...
# Comando para guppy_basecaller usando GPU
guppy_basecaller -r -i ${RAWDIR} -s ${BASECALLDIR} -c ${CONFIG} -x auto  --gpu_runners_per_device ${GPUPERDEVICE} --chunk_size ${CHUNCKSIZE} --chunks_per_runner ${CHUNKPERRUNNER} --verbose_logs
echo "Basecalling concluido com sucesso!"

# Comando 2 - Demultiplex & adapter removal (comum a todos workflows)
echo "Executando guppy_barcoder..."
mkdir -p "${DEMUXDIR}"
if [[ $WF -eq 2 ]];
then
	#bm5 WF2 sem headcrop 18 (uso cutadapt)
	guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends  --detect_mid_strand_barcodes --detect_mid_strand_adapter --trim_barcodes  
else
	#bm7 WF3 com headcrop 18
	guppy_barcoder -r -i "${BASECALLDIR}/pass" -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER}
fi

# Comando 3 - QC
#echo "Executando pycoQC..."
source activate ngs
pycoQC -q -f ${BASECALLDIR}/sequencing_summary.txt -b ${DEMUXDIR}/barcoding_summary.txt -o ${BASECALLDIR}/wf${WF}/${RUNNAME}_pycoqc.html --report_title ${RUNNAME} --min_pass_qual ${QSCORE}
# QC de arquivos .fastq
# fastqc <nome_do_arquivo.fastqc>

# Concatena todos arquivos .fastq de cada barcode em um arquivo .fastq único
[ ! -d "${DEMUXCATDIR}" ] && mkdir -vp ${DEMUXCATDIR}
for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
    [ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > ${DEMUXCATDIR}/${i}.fastq
done

# Comando 4 - cutadapt
echo "Executando cutadapt..."
[ ! -d "${CUTADAPTDIR}" ] && mkdir -vp ${CUTADAPTDIR}
for i in $(find "${DEMUXCATDIR}" -type f -exec basename {} .fastq \;); do
	cutadapt -g ${PRIMER} -e 0.2 --discard-untrimmed -o ${CUTADAPTDIR}/${i}.fastq ${DEMUXCATDIR}/${i}.fastq
done;

# Comando 5 - NanoFilt
echo "Executando NanoFilt..."
[ ! -d "${NANOFILTDIR}" ] && mkdir -vp ${NANOFILTDIR}
for i in $(find "${CUTADAPTDIR}" -type f -exec basename {} .fastq \;); do
	NanoFilt -l ${LENGTH} < ${CUTADAPTDIR}/${i}.fastq > ${NANOFILTDIR}/${i}.fastq 
done;

# Comando 6 - prinseq-lite.pl
# Link: https://chipster.csc.fi/manual/prinseq-complexity-filter.html
echo "Executando prinseq-lite.pl..."
[ ! -d "${PRINSEQDIR}" ] && mkdir -vp ${PRINSEQDIR}
for i in $(find "${NANOFILTDIR}" -type f -exec basename {} .fastq \;); do
	prinseq-lite.pl -fastq ${NANOFILTDIR}/${i}.fastq -out_good ${PRINSEQDIR}/${i}.good -out_bad ${PRINSEQDIR}/${i}.bad -graph_data ${PRINSEQDIR}/${i}.gd -no_qual_header -lc_method dust -lc_threshold 40
done;

# Converte arquivos .fastq em .fasta para query no blastn
[ ! -d "${QUERYDIR}" ] && mkdir -vp ${QUERYDIR}
for i in $(find "${PRINSEQDIR}"/*.good.fastq -type f -exec basename {} .good.fastq \;); do
	sed -n '1~4s/^@/>/p;2~4p' ${PRINSEQDIR}/${i}.good.fastq > ${QUERYDIR}/${i}.fasta
done;

# Comando 7 - nblast
# Script: makeblastdb_refseq.sh
# Concatena todas as REFSEQs num arquivo refseq.fasta único e cria o BLASTDB
# Script: fasta2acc.sh 
# Extrai do arquvio refseq.fasta a lista acesso refseq_acc.txt
# Script: acc2taxid.sh 
# Cria a partir do arquivo refser_acc.txt o arquivo refseq_map.txt que mapeia os taxid (números que identificam as espécies taxonômica)

# Busca as reads no BLASTDB local
echo "Executando blastn..."
[ ! -d "${BLASTDIR}" ] && mkdir -vp ${BLASTDIR}
for i in $(find "${QUERYDIR}" -type f -exec basename {} .fasta \;); do
	echo "Carregando os dados ${BLASTDIR}/${i}..."
	blastn -db ${BLASTDBDIR}/refseq -query ${QUERYDIR}/${i}.fasta -out ${BLASTDIR}/${i}.blastn -outfmt "6 sacc staxid" -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
	# Busca remota
	# blastn -db nt -remote -query ${QUERYDIR}/${i}.fasta -out ${BLASTDIR}/${i}.blastn -outfmt "6 qacc saccver pident sscinames length mismatch gapopen evalue bitscore"  -evalue 0.000001 -qcov_hsp_perc 90 -max_target_seqs 1
done;

# Comandos 8, 9  & 10 - minimap2, samtools & racon
echo "Executando minimap2, samtools & racon..."
[ ! -d "${READLEVELDIR}" ] && mkdir -vp ${READLEVELDIR}

# Cria o arquivo índice do genoma humano para reduzir o tempo de alinhamento
# minimap2 -d ${HUMANREFMMI} ${HUMANREFSEQ}

# Loop para analisar todos barcodes, um de cada vez
for i in $(find ${DEMUXCATDIR} -type f -name "*.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	# Alinha as reads contra o arquivo indice do genoma humano e ordena os segmentos
	minimap2 -ax map-ont -t ${THREADS} ${HUMANREFMMI} ${DEMUXCATDIR}/${i}.fastq | samtools sort -@ 12 -o ${READLEVELDIR}/${i}.sorted.bam -
	# Indexa o arquivo para acesso mais rápido
	samtools index -@ ${THREADS} ${READLEVELDIR}/${i}.sorted.bam
	# Filtra os segmento não mapeados Flag 4 (-f 4) para um novo arquivo .bam 
	samtools view -bS -f 4 ${READLEVELDIR}/${i}.sorted.bam > ${READLEVELDIR}/${i}.unmapped.sam -@ 12
	# Salva os dados no formato .fastq
	samtools fastq -f 4 ${READLEVELDIR}/${i}.unmapped.sam > ${READLEVELDIR}/${i}.unmapped.fastq -@ 12
	# Alinhar todas as reads com elas mesmas para produzir sequencias consenso a partir do overlap de reads
	minimap2 -ax ava-ont -t ${THREADS} ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.overlap.sam
	# Correção de erros a partir das sequencias consenso
	racon -t ${THREADS} -f -u ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.overlap.sam ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.corrected.fasta
done

fi # Fim do desvio para execução 

# Comando 11 - kraken2
echo "Executando o Kraken2..."
for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	echo "Carregando os dados ${READLEVELDIR}/${i}..."
	#kraken2 --db ${KRAKENDB} --threads ${THREADS} --report ${READLEVELDIR}/${i}_report.txt --report-minimizer-data --output ${READLEVELDIR}/${i}_output.txt ${READLEVELDIR}/${i}.corrected.fasta
	kraken2 --db ${KRAKENDB} --quick --threads ${THREADS} --report ${READLEVELDIR}/${i}_report.txt --output ${READLEVELDIR}/${i}_output.txt ${READLEVELDIR}/${i}.corrected.fasta
done

echo "Comando(s) executado(s) com sucesso!"
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
#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${CHIKVREFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.chikv.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.chikv.sorted.bam > ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${CHIKVREFSEQ} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.chikv -n N -i ${i}
#done

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
