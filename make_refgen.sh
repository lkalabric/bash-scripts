#!/bin/bash

# script: make_refgen.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# objetivo: extrai sequencias fasta do arquivo refseq.fasta usando a lista refgen e cria o arquivo refgen.fasta
# criação: 01 DEZ 2021
# ultima atualização: 01 DEZ 2021

REFSEQDIR="${HOME}/data/REFSEQ"

# Extrai os genomas referencia do arquivo refseq.fasta a partir da lista disponível em refgen
# Link: https://www.biostars.org/p/319099/
cat ${REFSEQDIR}/refgen | cut -f 2 | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${REFSEQDIR}/refseq.fasta > ${REFSEQDIR}/refgen.fasta
# grep -A 1 -wFf $(cut ${REFSEQDIR}/refgen -f 2) ${REFSEQDIR}/refseq.fasta > ${REFSEQDIR}/refgen.fasta

