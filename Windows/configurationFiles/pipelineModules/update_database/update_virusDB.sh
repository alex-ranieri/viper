#!/bin/bash

# Updating pangolin
echo "Updating pangolin data"
pangolin --update-data

# Update SARS-CoV-2 Nextclade data
echo "Updating SARS-CoV-2 Nextclade data"
nextclade dataset get --name='sars-cov-2' --output-dir="$PIPELINE/SARS-CoV-2/nextstrain_files/"

# Updating Influenza Nextclade data

echo "Updating Influenza Nextclade data"
nextclade dataset get --name='flu_h1n1pdm_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/H1"
nextclade dataset get --name='flu_h3n2_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/H3"
nextclade dataset get --name='flu_vic_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/Vic"
nextclade dataset get --name='flu_yam_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/Yam"

# Update DENV data
echo "Updating Dengue Nextclade data"
REPO="alex-ranieri/denvLineages"
FOLDERS=("DENV1" "DENV2" "DENV3" "DENV4")
DOWNLOAD_PATH="$PIPELINE/DENV/nextclade_files"

# Function to check if a directory exists and delete it if it does
delete_directory() {
    if [ -d "$1" ]; then
        echo "Deleting existing directory: $1"
        rm -rf "$1"
    fi
}

# Fetch the GitHub API response for each folder and download its contents
for FOLDER in "${FOLDERS[@]}"; do
    # Delete existing directory if it exists
    delete_directory "$DOWNLOAD_PATH/denv${FOLDER: -1}"

    # Create directory for the current folder
    mkdir -p "$DOWNLOAD_PATH/denv${FOLDER: -1}"

    # Fetch the GitHub API response for the folder
    API_URL="https://api.github.com/repos/$REPO/contents/Nextclade_V3_data/$FOLDER"
    FILES=$(curl -s "$API_URL" | grep -Eo '"download_url": "[^"]+"' | cut -d '"' -f 4)

    # Loop through each file and download it
    for FILE in $FILES; do
        FILENAME=$(basename "$FILE")
        echo "Downloading $FILENAME into $DOWNLOAD_PATH/denv${FOLDER: -1}..."
        curl -s -o "$DOWNLOAD_PATH/denv${FOLDER: -1}/$FILENAME" -J -L "$FILE"
    done
done

echo "Done"