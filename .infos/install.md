
## list of some dependencies 

## some are still missing here!

**Protint** 


direct install: 

```sh
wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz 
tar zxvf ProtHint-2.6.0.tar.gz
cd ProtHint-2.6.0/bin 
then add to ~/.bashrc
echo -e "\n#Path to ProtHint\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../../
```

**Diamond**

direct install: 

```sh
mkdir diamond ; cd diamond
wget https://github.com/bbuchfink/diamond/releases/download/v2.1.1/diamond-linux64.tar.gz`
tar zxvf diamond-linux64.tar.gz

then add to ~/.bashrc
echo -e "\n#Path to diamond\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../

```

**cdbfasta**:

direct install: 

```
git clone https://github.com/gpertea/cdbfasta.git`
cd cdbfasta
make all 
then add to ~/.bashrc
echo -e "\n#Path to cdbfasta\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  

```


**Augustus**

see details at https://github.com/Gaius-Augustus/Augustus
with a bit of luck this could work:

direct install: 

```
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus
make augustus
```

For augustus I also had to install bamtools: 
bamtools : https://github.com/pezmaster31/bamtools
some change may be necessary to install bamtools properly....

** boost ** was also necessary on some cluster ....


the rest was OK but you may need to carrefully read the Readme of Augustus to make it work


### add augsutus script/config/bin to the .bashrc ###

```
export AUGUSTUS_CONFIG_PATH=/home/path/to/augustus/config
export AUGUSTUS_BIN_PATH=/home/path/to/augustus/bin/
export AUGUSTUS_SCRIPTS_PATH=/home/path/to/augustus_scripts
```


**genemark** 

wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_pxuuc/gmes_linux_64.tar.gz

the gm_key is necessary
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_pxuuc/gm_key_64.gz

it must be copied to your home
```sh
gunzip gm_key_64.gz 
mv gm_key_64 ~/.gm_key
```

Note: you must register [online](http://exon.gatech.edu/GeneMark/license_download.cgi) for genemark


#Other dependencies: 

**TSEBRA** available [here](https://github.com/Gaius-Augustus/TSEBRA)

direct install: 

```sh
git clone https://github.com/Gaius-Augustus/TSEBRA 
cd TSEBRA/bin
path=$(pwd)
echo -e "\n#Path to TSEBRA\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../
```


**trimmomatic** software available [here](http://www.usadellab.org/cms/?page=trimmomatic)

direct install:
```sh
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
```


**GSNAP** for read mappig available [here](http://research-pub.gene.com/gmap/)

direct install :
```sh
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

```


[repeatmodeler](https://www.repeatmasker.org/RepeatModeler/) and [repeatmasker](https://www.repeatmasker.org/)


**gffread** to extract CDS from fasta [see](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)

direct install: 

```sh
git clone https://github.com/gpertea/gffread
cd gffread
make release
echo -e "\n#Path to gffread\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../
```  

**transeq**  from EMBOSS to convert fasta into protein [click_here](https://www.bioinformatics.nl/cgi-bin/emboss/help/transeq) 

direct install: 

```sh
#to install we advocate to use mamba - an equivlent to conda but faster:
mambe create -p braker_env emboss=6.6.0
```

note: the same results could be obtained with gffread instead

**BUSCO** for quality assesment (https://busco.ezlab.org/)

direct install: 

```sh
mamba activate braker_env 
mamba install -c bioconda busco=5.5.0
```

**orthofinder** : 

direct install: 

```sh
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

```

**MCScanX** : 

direct install: 

```sh 
git clone https://github.com/wyp1125/MCScanX
cd MSCanX ; make -j 8
#if command was successfull then add to path:
path=$(pwd)
#add it to the bashrc
echo -e "\n#Path to MScanX\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
```

**muscle** 

direct install: 

```sh
mkdir muscle ; cd muscle
wget wget https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.linux_intel64
ln -s muscle5.1.linux_intel64 muscle
chmod +x muscle
path=$(pwd)
echo -e "\n#Path to muscle\n export PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../
```

**yn00/paml:** 


direct install: 

```sh
wget https://github.com/abacus-gene/paml/releases/download/4.10.7/paml-4.10.7-linux-X86_64.tgz
tar zxvf paml-4.10.7-linux-X86_64.tgz
cd paml-4.10.7/bin
path=$(pwd)
echo -e "\n#Path to paml\n export PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../../
```
    
**translatorx_vLocal.pl** 

direct install: 

```sh
wget http://161.111.160.230/cgi-bin/translatorx_vLocal.pl
path=$(pwd)
echo -e "\n#Path to translatorx\n export PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
```

**minimap2** 

direct install: 

```sh
curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
cd ./minimap2-2.26_x64-linux 
path=$(pwd)
echo -e "\n#Path to minimap2\n export PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
cd ../
```


**samtools**

direct install: 

```sh
mamba activate braker_env 
mamba install -c bioconda samtools samtools=1.18
```

**R** with several packages 


**optional:** (not used here)

**miniprot**

**interproscan** for [annotation](https://interproscan-docs.readthedocs.io/en/latest/index.html):  

**[agat](https://agat.readthedocs.io/en/latest/index.html)** for statistics, quality assesment, errors

