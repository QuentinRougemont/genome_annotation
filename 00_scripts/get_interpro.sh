#!/bin/bash
#QR
#2024

#get interpro-scan
    
mkdir my_interproscan 2>/dev/null
cd my_interproscan || exit
wget -q https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/interproscan-5.71-102.0-64-bit.tar.gz
wget -q https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/interproscan-5.71-102.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
#md5sum -c interproscan-5.71-102.0-64-bit.tar.gz.md5

tar -pxvzf interproscan-5.71-102.0-*-bit.tar.gz
cd interproscan-5.71-102.0 || exit

python3 setup.py -f interproscan.properties
    
path=$(pwd)
echo -e "\n#Path to interproscan\nexport PATH=\$PATH:$path" >> ~/.bashrc
source ~/.bashrc

#maybe set interpro to false in the config file if installation fails 
