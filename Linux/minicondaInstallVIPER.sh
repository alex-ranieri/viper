#!/bin/bash

# Captura o nome do usuário que está executando o script
usuario=$USER

# Comando original com substituição da variável
comando="wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && bash /tmp/miniconda.sh -b -p /home/$usuario/miniconda3 && rm /tmp/miniconda.sh"

# Executa o comando
eval $comando

# Adiciona o Miniconda ao PATH do usuário
echo "export PATH=/home/$usuario/miniconda3/bin:\$PATH" >> /home/$usuario/.bashrc

# Configura o Conda para inicialização automática
echo "source /home/$usuario/miniconda3/etc/profile.d/conda.sh" >> /home/$usuario/.bashrc

# Inicializa o Conda automaticamente
/home/$usuario/miniconda3/bin/conda init bash
source /home/$usuario/.bashrc

# Configurações do ambiente Conda para o usuário
echo "Configuring Conda env to user"
/home/$usuario/miniconda3/bin/conda config --add channels conda-forge
/home/$usuario/miniconda3/bin/conda config --add channels anaconda
/home/$usuario/miniconda3/bin/conda config --add channels bioconda
/home/$usuario/miniconda3/bin/conda config --set channel_priority disabled

# Criação do ambiente Conda "GenomeAssembler" para o novo usuário
echo "Creating Conda env 'VIPER'"
/home/$usuario/miniconda3/bin/conda install -y -c conda-forge micromamba
echo "Installing packages in 'VIPER' env using micromamba"
/home/$usuario/miniconda3/bin/micromamba create -y -f VIPER.yml
