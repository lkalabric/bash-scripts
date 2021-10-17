# Clone remote to local
#git clone https://github.com/lkalabric/bash-scripts.git /home/brazil1/bash-scripts

# Define nome de usuário do Git
cd ~/bash-scripts
git config --global user.name "lkalabric"
git config --global user.email "luciano.kalabric@fiocruz.br"

# Configura o método de autenticação do usuário
gh user lkalabric

# Commit local
git commit -a -m "Atualização do script assembly.sh."

# Push changes to the server
git push
