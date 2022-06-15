#!/bin/bash

# Cria o banco de dados BLAST a partir de arquivos .fasta
# Autor: Luciano Kalabric Silva
# Data: 28/08/2021
# Referencia: https://www.ncbi.nlm.nih.gov/books/NBK569841/

REFSEQDIR=${HOME}/data/REFSEQ
BLASTDBDIR=${HOME}/data/BLAST_DB

[[ -r ${REFSEQDIR}/refseq.fasta ]] && rm ${REFSEQDIR}/refseq.fasta
[[ -d ${BLASTDBDIR} ]] && rm -r ${BLASTDBDIR}
[ ! -d ${BLASTDBDIR} ] && mkdir -vp ${BLASTDBDIR}

# Concatena todos os arquivos .fasta em refseq.fasta
echo "Concatenando as sequencias referência em refseq.fasta..."
if [[ ! -f ${REFSEQDIR}/refseq.fasta ]]
then
	find ${REFSEQDIR} -type f -name '*.fasta' -print0 | sort -z | xargs -0 cat > "${REFSEQDIR}/refseq.fasta"
	# mv refseq.fasta ${REFSEQDIR}
fi

# Processando a linha de descrição para conter apenas o número de acesso sem espaços
echo "Processando os labels do arquivo refseq.fasta..."
cp ${REFSEQDIR}/refseq.fasta ${REFSEQDIR}/refseq.old
while read -r line; do
	if echo "$line" | grep ">";
	then
    	echo "$line" | cut -d "." -f 1 >>${REFSEQDIR}/refseq.temp
	else
		echo "$line" >>${REFSEQDIR}/refseq.temp
	fi
done < ${REFSEQDIR}/refseq.fasta
mv ${REFSEQDIR}/refseq.temp ${REFSEQDIR}/refseq.fasta

# Cria a lista de números de acc Genbank a partir do arquivo .fasta
echo "Criando o arquivo refseq.acc..."
[[ -f ${REFSEQDIR}/refseq.acc ]] && rm  ${REFSEQDIR}/refseq.acc
grep ">" ${REFSEQDIR}/refseq.fasta | sed 's/>//' | cut -d " " -f 1 > ${REFSEQDIR}/refseq.acc

# Cria a lista de taxid a partir nos números de acc Genbank
[[ -f ${REFSEQDIR}/refseq.map ]] && rm  ${REFSEQDIR}/refseq.map
# Retrive Taxid
echo "Criando o arquivo refseq.map..."
while read -r line; do
	echo "$line "$(efetch -db nuccore -id "$line" -format docsum | xtract -pattern DocumentSummary -element TaxId) >>${REFSEQDIR}/refseq.map
done < ${REFSEQDIR}/refseq.acc

# Cria o banco de dados Blast a partir do arquivo refseq.fasta
echo "Criando o banco de dados BLAST_DB/refseq..."
makeblastdb -in ${REFSEQDIR}/refseq.fasta -parse_seqids -blastdb_version 5 -taxid_map ${REFSEQDIR}/refseq.map -dbtype nucl -out ${BLASTDBDIR}/refseq
echo "Banco de dados criado com sucesso!"

# Faz o donwload do taxdb
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz -P ${BLASTDBDIR}
