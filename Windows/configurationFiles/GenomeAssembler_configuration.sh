#!/bin/bash

# Configuração dos servidores DNS públicos do Google
echo "nameserver 8.8.8.8" >> /etc/resolv.conf
echo "nameserver 8.8.4.4" >> /etc/resolv.conf

# Atualiza a lista de pacotes 
echo "Updating list of packages..."
apt-get update

# Instala pacotes essenciais
echo "Installing essential packages..."
apt install -y build-essential wget

# Instala o neofetch
apt install -y neofetch bc csvkit libgl1

# Obtém o nome do usuário como argumento
password=$1

# Cria o novo usuário e solicita uma senha
echo "Creating VIPER user"
#useradd -m -s /bin/bash viper
#echo -e "$password\n$password\n" | passwd viper
useradd -m -s /bin/bash -p $(echo "$password" | openssl passwd -1 -stdin) viper > /dev/null 2>&1

# Define o novo usuário como usuário padrão em WSL
echo "root:viper" | chpasswd
echo "viper ALL=(ALL) ALL" >> /etc/sudoers
echo "export WSL_INTEROP=viper" >> /etc/environment
echo "[user]" >> /etc/wsl.conf
echo "default=viper" >> /etc/wsl.conf

# Ajusta os glitchs de visualização de interface gráfica
echo "[system-distro-env]" >> /etc/wsl.conf
echo ";disable GPU in system-distro" >> /etc/wsl.conf
echo "LIBGL_ALWAYS_SOFTWARE=1" >> /etc/wsl.conf

# Baixa e instala o Miniconda3 na pasta home do novo usuário
echo "Downloading and installing Miniconda3"
sudo -u viper bash -c "wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && bash /tmp/miniconda.sh -b -p /home/viper/miniconda3 && rm /tmp/miniconda.sh"

# Adiciona o Miniconda ao PATH do novo usuário
echo "export PATH=/home/viper/miniconda3/bin:\$PATH" >> /home/viper/.bashrc

# Configura o Conda para inicialização automática
echo "source /home/viper/miniconda3/etc/profile.d/conda.sh" >> /home/viper/.bashrc

# Inicializa o Conda automaticamente
sudo -u viper bash -c "/home/viper/miniconda3/bin/conda init bash"
sudo -u viper bash -c "source /home/viper/.bashrc"

# Configurações do ambiente Conda para o novo usuário
echo "Configuring Conda env to new user"
sudo -u viper bash -c "/home/viper/miniconda3/bin/conda config --add channels conda-forge"
sudo -u viper bash -c "/home/viper/miniconda3/bin/conda config --add channels anaconda"
sudo -u viper bash -c "/home/viper/miniconda3/bin/conda config --add channels bioconda"
sudo -u viper bash -c "/home/viper/miniconda3/bin/conda config --set channel_priority disabled"

# Criação do ambiente Conda "GenomeAssembler" para o novo usuário
echo "Creating Conda env 'VIPERGenomeAssembler'"
sudo -u viper bash -c "/home/viper/miniconda3/bin/conda install -y -c conda-forge micromamba" 
echo "Installing packages in Conda env 'VIPERGenomeAssembler' using micromamba"
sudo -u viper bash -c "/home/viper/miniconda3/bin/micromamba create -y --prefix /home/viper/miniconda3/envs/VIPERGenomeAssembler python=3.8 fastqc trimmomatic fastp cutadapt spades bowtie2 samtools ivar pangolin minimap blast quast pilon nextclade seqtk sra-tools bandage iqtree xlrd openpyxl pyqt"

# Exportando comando para ativar o ambiente GenomeAssembler por padrão para o usuário
echo "conda activate VIPERGenomeAssembler" >> /home/viper/.bashrc 

# Exporta neofetch
echo "neofetch" >> /home/viper/.bashrc

# Mensagem de conclusão
echo "User viper created successfully, Miniconda3 installed in the home folder, Conda environment configured and packages installed."

# Configure VIPER modules PATH
echo "Configuring pipeline modules in /home/viper/pipelineModules"
sudo -u viper bash -c "cp -r pipelineModules /home/viper/"
sudo -u viper bash -c "cp -r viperGUI /home/viper/"

# Adiciona o diretório ao PATH
PIPELINE_PATH=/home/viper/pipelineModules
echo "export PATH=\$PATH:$PIPELINE_PATH" >> /home/viper/.bashrc
# Define a variável PIPELINE
echo "export PIPELINE=$PIPELINE_PATH" >> /home/viper/.bashrc
ln -s $PIPELINE_PATH/SARS-CoV-2/Exec_COV_assembly_pipeline_Illumina_CeVIVAS.sh /usr/local/bin/VIPER_CoV.sh
ln -s $PIPELINE_PATH/DENV/Exec_DENV_assembly_pipeline_Illumina.sh /usr/local/bin/VIPER_DENV.sh
ln -s $PIPELINE_PATH/Influenza/Exec_FLU_assembly_pipeline_Illumina.sh /usr/local/bin/VIPER_Influenza.sh
ln -s /home/viper/viperGUI/RunVIPER.bash /usr/local/bin/RunVIPER.bash


sudo -u viper bash -c "bash -l"

echo "Configuration complete! You can now close this window!"