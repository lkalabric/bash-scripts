#!/bin/bash

# script: metagenomics_wf1.sh
# autores: Laise de Moraes <laisepaixao@live.com> & Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# criação: 25 AGO 2021
# última atualização: 25 NOV 2021
# versão wf1: definição de variáveis, criação dos diretórios de I/O e basecalling

# Validação da entrada de dados na linha de comando
RUNNAME=$1 	# Nome do dado passado na linha de comando
MODEL=$2	# Modelo de basecalling fast hac sup
WF=1		# Workflow de bioinformatica 1, 2 ou 31

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

# Parâmetros Guppy basecaller (ONT)
CONFIG="dna_r9.4.1_450bps_${MODEL}.cfg" #dna_r9.4.1_450bps_fast.cfg dna_r9.4.1_450bps_hac.cfg dna_r9.4.1_450bps_sup.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
 
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

summary1 () {
# Step 0 - Sumario do sequenciamento (dados disponíveis no arquivo report*.pdf)
echo "Sumário da corrida"
echo "Total files:"
ls $(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;) | wc -l
echo "Total reads:"
# h5ls "$(find ${RAWDIR} -type f -name "*.fast5" -exec dirname {} \;)"/*.fast5 | wc -l
}
# Estima o tempo da execução
export -f summary1
echo 'summary1' | /usr/bin/time -o "~/performance-analysis/${RUNNAME}_summary1.time" /bin/bash

exit 0

# Pausa a execução para debug
# read -p "Press [Enter] key to continue..."

basecalling () {
# Step 1 - Basecalling (comum a todos workflows)
# Esta etapa pode ser realizada pelo script guppy_gpu_v1_ag.sh no LAPTOP-Yale
if [ ! -d $BASECALLDIR ]; then
	echo -e "\nExecutando guppy_basecaller..."
	mkdir -vp $BASECALLDIR
	# Esta etapa está sendo realizada pelo script guppy_gpu_v1_ag.sh no LAPTOP-Yale
	# Comando para guppy_basecaller usando GPU
	guppy_basecaller -r -i ${RAWDIR} -s "${BASECALLDIR}" -c ${CONFIG} -x auto  --gpu_runners_per_device ${GPUPERDEVICE} --chunk_size ${CHUNCKSIZE} --chunks_per_runner ${CHUNKPERRUNNER} --verbose_logs
fi
}
# Estima o tempo da execução
export -f basecalling
echo 'basecalling' | /usr/bin/time -o "~/performance-analysis/${RUNNAME}_basecalling.time" /bin/bash
