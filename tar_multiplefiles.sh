#!/bin/bash

# script: tar_multiplefiles.sh
# autor: Luciano Kalabric <luciano.kalabric@fiocruz.br>
# instituição: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# objetivo: compacta multiplos arquivos de um diretório em multiplos arquivos .tar
# criação: 08 JUL 2022
# ultima atualização: 08 JUL 2022
# atualização: criação do arquivo

# solução 1
for file in `ls $1`; do tar -czvf $file.tar.gz $file ; done

# solução 2
# find . -type f -execdir zip '{}.zip' '{}' \;
