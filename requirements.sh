

#list all requirements for installation

#install mamba for linux if not already installed!

command='mamba'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "please install mamba first or you'll get trouble"
    exit 1
fi

mamba create -p braker_env  -c anaconda perl biopython
mamba install -c bioconda perl-app-cpanminus perl-hash-merge perl-parallel-forkmanager \
    perl-scalar-util-numeric perl-yaml perl-class-data-inheritable \
    perl-exception-class perl-test-pod perl-file-which  perl-mce \
    perl-threaded perl-list-util perl-math-utils cdbtools \
    perl-list-moreutils

mamba activate braker_env
mamba install -c bioconda perl-file-homedir


#**Protint** 
#direct install: 

wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz 
tar zxvf ProtHint-2.6.0.tar.gz
cd ProtHint-2.6.0/bin 
then add to ~/.bashrc
echo -e "\n#Path to ProtHint\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../../


#**Diamond**
#direct install: 

mkdir diamond ; cd diamond
wget https://github.com/bbuchfink/diamond/releases/download/v2.1.1/diamond-linux64.tar.gz`
tar zxvf diamond-linux64.tar.gz

then add to ~/.bashrc
echo -e "\n#Path to diamond\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../


#**cdbfasta**:
#direct install: 
git clone https://github.com/gpertea/cdbfasta.git`
cd cdbfasta
make all 
then add to ~/.bashrc
echo -e "\n#Path to cdbfasta\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  


#**bamtools**
bamtools : https://github.com/pezmaster31/bamtools
git clone https://github.com/pezmaster31/bamtools
cd bamtools
mkdir build
cd build
path=$(pwd)
cmake -DCMAKE_INSTALL_PREFIX=$(pwd) ..

make -j8 

make DESTDIR=$path install
cd ../build/src
path=$(pwd)
echo -e "\n#Path to bamtools\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  


#**Augustus**
#direct install: 
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus
make augustus
```
### add augsutus script/config/bin to the .bashrc ###
#note you may encounter several error and have to comment/uncomment or change path in "common.mk" especially without root privilege
cd 
export AUGUSTUS_CONFIG_PATH=/home/path/to/augustus/config
export AUGUSTUS_BIN_PATH=/home/path/to/augustus/bin/
export AUGUSTUS_SCRIPTS_PATH=/home/path/to/augustus_scripts
```


#**genemark** 
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_pxuuc/gmes_linux_64.tar.gz
echo "to get genemark to work you must register online at "http://exon.gatech.edu/GeneMark/license_download.cgi"
echo "the gm_key is necessary" 
echo "please download the the GeneMark-ES/ET/EP under the appropriate linux kernel"
echo "it must be copied to your home
echo -e "do the following\n
gunzip gm_key_64.gz 
mv gm_key_64 ~/.gm_key "

#**TSEBRA** available [here](https://github.com/Gaius-Augustus/TSEBRA)
#direct install: 
git clone https://github.com/Gaius-Augustus/TSEBRA 
cd TSEBRA/bin
path=$(pwd)
echo -e "\n#Path to TSEBRA\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../


#**trimmomatic** software available [here](http://www.usadellab.org/cms/?page=trimmomatic)
#direct install:
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip


#**GSNAP** for read mappig available [here](http://research-pub.gene.com/gmap/)
#direct install :
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2023-10-10.v2.tar.gz
tar zxf gmap-gsnap-2023-10-10.v2.tar.gz
cd gmap-2023-10-10/
path=$(pwd)
./configure --prefix=$path  
make -j8
make install

cd bin
echo -e "\n#Path to gsnap-gmap\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  


#[repeatmodeler](https://www.repeatmasker.org/RepeatModeler/) and [repeatmasker](https://www.repeatmasker.org/)
#to do!


#**gffread** to extract CDS from fasta [see](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
#direct install: 
git clone https://github.com/gpertea/gffread
cd gffread
make release
echo -e "\n#Path to gffread\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../

#**transeq**  from EMBOSS to convert fasta into protein [click_here](https://www.bioinformatics.nl/cgi-bin/emboss/help/transeq) 
#direct install: 
mamba activate braker_env 
mambe create -p braker_env emboss=6.6.0


#**BUSCO** for quality assesment (https://busco.ezlab.org/)
#direct install: 
mamba activate braker_env 
mamba install -c bioconda busco=5.5.0

#**orthofinder** : 
#direct install: 
wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder.tar.gz
tar zxf OrthoFinder.tar.gz
cd OrthoFinder
#if command was successfull then add to path:
path=$(pwd)
echo -e "\n#Path to othofinder\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
    
#also add diamond, fastme and mcl which are present within orthofinder:
cd bin/
path=$(pwd)
echo -e "\n#Path to diamond\n export PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../../


#**MCScanX** : 
#direct install: 
git clone https://github.com/wyp1125/MCScanX
cd MSCanX ; make -j 8
#if command was successfull then add to path:
path=$(pwd)
#add it to the bashrc
echo -e "\n#Path to MScanX\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  

#**muscle** 
#direct install: 
mkdir muscle ; cd muscle
wget wget https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.linux_intel64
ln -s muscle5.1.linux_intel64 muscle
chmod +x muscle
path=$(pwd)
echo -e "\n#Path to muscle\n export PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../

#**yn00/paml:** 
#direct install: 

wget https://github.com/abacus-gene/paml/releases/download/4.10.7/paml-4.10.7-linux-X86_64.tgz
tar zxvf paml-4.10.7-linux-X86_64.tgz
cd paml-4.10.7/bin
path=$(pwd)
echo -e "\n#Path to paml\n export PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../../
    
#**translatorx_vLocal.pl** 
#direct install: 
wget http://161.111.160.230/cgi-bin/translatorx_vLocal.pl
path=$(pwd)
echo -e "\n#Path to translatorx\n export PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  

#**minimap2** 
#direct install: 
curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
cd ./minimap2-2.26_x64-linux 
path=$(pwd)
echo -e "\n#Path to minimap2\n export PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../


#**samtools**
#direct install: 
mamba activate braker_env 
mamba install -c bioconda samtools samtools=1.18

#**R** with several packages 
mamba activate braker_env 
mamba install -c conda-forge r-base
