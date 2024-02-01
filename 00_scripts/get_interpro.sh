#QR
#2024

#get interpro-scan
    
mkdir my_interproscan
cd my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.66-98.0/interproscan-5.66-98.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.66-98.0/interproscan-5.66-98.0-64-bit.tar.gz.md5
    
#md5sum -c interproscan-5.66-98.0-64-bit.tar.gz.md5
 
tar -pxvzf interproscan-5.66-98.0-*-bit.tar.gz
cd interproscan-5.66-98.0

python3 setup.py -f interproscan.properties
    
path=$(pwd)
echo -e "\n#Path to interproscan\nexport PATH=\$PATH:$path" >> ~/.bashrc
source ~/.bashrc

#maybe set interpro to false in the config file if installation fails 
