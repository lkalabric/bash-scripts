#!/bin/bash

# Cria o banco de dados GMAP a partir de arquivos .fasta (.fa.gz)
# Autor: Luciano Kalabric Silva
# Data: 11/10/2022
# Última atualização: 11/10/2022
# Referencia: https://github.com/juliangehring/GMAP-GSNAP/blob/master/README

# Validação da entrada de dados na linha de comando
GENOMA=$1 	# Nome do genoma
FASTA=$2      	# Arquivo FASTA
FLAGS=$3	# Flags, se aplicável
if [[ $# -eq 0 ]]; then
	echo "Falta o nome do genoma referência ou do arquivo(s) FASTA(s)!"
	echo "Sintáxe: ./make_gmapdb.sh <GENOMA> <FASTA>"
  	echo "Exemplo: ./make_gmapdb.sh GRCh38 data/GRCh38/GRCh38.p13.genome.fa.gz -g"
  	echo "Exemplo: ./make_gmapdb.sh HCV data/REFGEN/hcv1a.fasta"
	exit 0
fi

# Diretório do banco de dados 
GMAPDBDIR=${HOME}/data/GMAPDB

# Reseta o diretório antes de criar um novo banco de dados
[[ -d ${GMAPDBDIR}/${GENOMA} ]] && rm -r ${GMAPDBDIR}/${GENOMA}

# Monta o banco de dados
gmap_build -d ${GENOMA} ${FASTA} ${FLAGS}
