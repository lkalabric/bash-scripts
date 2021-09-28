# O arquivo taxin contém a lista de taxons de interesse
# Atualizar na medida do necessário
RUNNAME=$1
WF=31
BARCODE=$2
RESULTSDIR="$HOME/ngs-analysis/$RUNNAME/wf$WF"
READSLEVELDIR="${RESULTSDIR}/READS_LEVEL"

if [[ $# -ne 2 ]]; then
	echo "Faltando algum parâmetro!"
	echo "Sintáxe: ./kraken2_quick_report.sh <BIBLIOTECA_MODEL> <BARCODE>"
	echo "Exemplo: ./kraken2_quick_report.sh DENV_FTA_1_hac barcode01"
	exit 0
fi
while read -r line ; do
	count=$(agrep -q -w "$line" ${READSLEVELDIR}/${BARCODE}_report.txt | cut -f 2)
	echo "$line - $count"
done < <(cat /home/brazil1/data/REFSEQ/taxin | tr '\t' ';')
