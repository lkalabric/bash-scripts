#!/bin/bash

# Cria o banco de dados GMAP a partir de arquivos .fasta (.fa.gz)
# Autor: Luciano Kalabric Silva
# Data: 11/10/2022
# Última atualização: 11/10/2022
# Referencia: https://github.com/juliangehring/GMAP-GSNAP/blob/master/README

# Validação da entrada de dados na linha de comando
GENOMA=$1 	  # Nome do genoma
FASTA=$2      # Arquivo FASTA
if [[ $# -eq 0 ]]; then
	echo "Falta o nome do genoma referência ou do arquivo(s) FASTA(s)!"
	echo "Sintáxe: ./make_gmapdb.sh <GENOMA> <FASTA>"
  echo "Exemplo: ./make_gmapdb.sh GRCh38 GRCh38.p13.genome.fa.gz
	exit 0
fi

# Diretórios de dados
FASTADIR=${HOME}/data/GENOMA/${1}

# Diretório do banco de dados 
GMAPDBDIR=${HOME}/data/GMAPDB

# Reseta o diretório antes de criar um novo banco de dados
[[ -d ${GMAPDBDIR}/${1} ]] && rm -r ${GMAPDBDIR}/${1}

# Monta o banco de dados
gmap_build -d $1 ${FASTADIR}/GRCh38.p13.genome.fa.gz -g
