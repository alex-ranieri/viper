#!/bin/bash

# Obtém o caminho absoluto do diretório atual
PIPELINE_PATH=$(realpath .)

# Adiciona o diretório ao PATH
echo "export PATH=\$PATH:$PIPELINE_PATH" >> ~/.bashrc

# Carrega as alterações no shell atual
source ~/.bashrc

echo "Current directory added to PATH as PIPELINE"