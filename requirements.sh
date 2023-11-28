

#list all requirements for installation - attempts to download and install missing deps
#very minimal requirement than will not be installed here:
#wget
#curl
#git 
#make 
#cmake
#make 
#mamba 
#java

command='wget'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi

command='curl'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi

command='make'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi

command='cmake'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    echo "this can be done with mamba - see below"
    echo "once mamba is installed use: mamba install -c conda-forge cmake" 
    exit 1
fi


#install mamba for linux if not already installed!
command='mamba'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "please install mamba first or you'll get trouble"
    exit 1
fi

mamba create -p braker_env  -c anaconda perl biopython
path=$(pwd)
eval "$(conda shell.bash hook)"
conda activate $path/braker_env
mamba install -c bioconda perl-app-cpanminus perl-hash-merge perl-parallel-forkmanager \
    perl-scalar-util-numeric perl-yaml perl-class-data-inheritable \
    perl-exception-class perl-test-pod perl-file-which  perl-mce \
    perl-threaded perl-list-util perl-math-utils cdbtools \
    perl-list-moreutils
mamba install -c bioconda perl-file-homedir


#----------- a few tools to install directly with mamba :---------------------- 
#**samtools**
#direct install: 
command='samtools'
if ! command -v $command &> /dev/null
then
   #direct install: 
   mamba install -c bioconda samtools samtools=1.18
fi 

#**R** with several packages 
command='R'
if ! command -v $command &> /dev/null
then
   #direct install: 
   mamba install -c conda-forge r-base
fi 

#**transeq**  from EMBOSS to convert fasta into protein [click_here](https://www.bioinformatics.nl/cgi-bin/emboss/help/transeq) 
#direct install: 
command='transeq'
if ! command -v $command &> /dev/null
then
   #direct install: 
   mambe create -p braker_env emboss=6.6.0
fi
 
conda deactivate

#---------------------- OTHER TOOLS ---------------------------------------#
# test each command one by one and install them if necessary:
mkdir softs
cd softs 

#braker: 
command='braker.pl'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation" 
    git clone https://github.com/Gaius-Augustus/BRAKER
    cd BRAKER/scripts 
    chmod a+x *.pl *.py
    path=$(pwd)
    echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../../
fi

#**Protint** 
#direct install: 
command='ProtHint'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz 
    tar zxvf ProtHint-2.6.0.tar.gz
    cd ProtHint-2.6.0/bin 
    #then add to ~/.bashrc
    path=$(pwd)
    echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../../
fi


#**Diamond**
#direct install: 
command='diamond'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    mkdir diamond ; cd diamond
    wget https://github.com/bbuchfink/diamond/releases/download/v2.1.1/diamond-linux64.tar.gz
    tar zxvf diamond-linux64.tar.gz
    #then add to ~/.bashrc
    path=$(pwd)
    echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../ 

fi


#**cdbfasta**:
command='cdbfasta'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    #direct install: 
    git clone https://github.com/gpertea/cdbfasta.git
    cd cdbfasta
    make all 
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
	cd ../
    else
        echo -e "\n#ERROR : Installation failed please check everything"  
	exit 1
    fi
fi


#**bamtools**
command='bamtools'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    #direct install: 
    bamtools : https://github.com/pezmaster31/bamtools
    git clone https://github.com/pezmaster31/bamtools
    cd bamtools
    mkdir build
    cd build
    path=$(pwd)
    cmake -DCMAKE_INSTALL_PREFIX=$(pwd) ..
    make -j8 
    make DESTDIR=$path install
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
	cd ../build/src
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../../../
    else
       echo installation failed\nmake sur to have git, make and cmake
       exit 1
    fi
fi


#**Augustus**
command='Augustus'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will attempt a manual install" 
    echo "you may encounter several error and have to comment/uncomment or change path in "common.mk" especially without root privilege"
    git clone https://github.com/Gaius-Augustus/Augustus.git
    cd Augustus
    sed -i 's/COMPGENEPRED = true/COMPGENEPRED = false/g' common.mk
    sed -i 's/ZIPINPUT = true/ZIPINPUT = false\nBOOST = false/g' common.mk
    bt=$(command -v bamtools)
    echo "INCLUDE_PATH_BAMTOOLS     := -I$bt/src" >> common.mk
    echo "LIBRARY_PATH_BAMTOOLS     := -L$bt/usr/local/lib -Wl,-rpath,$bt/usr/local/lib" >> common.mk
    #attempt to make augustus:
    make augustus
    #if all was succesffull: ""
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
        path=$(pwd)
	echo -e "#path to AUGUSTUS:" >> ~/.bashrc 
        echo -e "export AUGUSTUS_CONFIG_PATH=$path/config " >> ~/.bashrc
        echo -e "export AUGUSTUS_BIN_PATH=$path/bin/ " >> ~/.bashrc
        echo -e "export AUGUSTUS_SCRIPTS_PATH=/$path/augustus_scripts " >> ~/.bashrc
	source ~/.bashrc
        cd ../
    else
	    echo "installation of augustus failed, look at the logs ! "
	    echo "see details here: https://github.com/Gaius-Augustus/Augustus/blob/master/docs/INSTALL.md"
	    exit 1
    fi
fi


#**genemark** 
#wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_pxuuc/gmes_linux_64.tar.gz
echo "to get genemark to work you must register online at http://exon.gatech.edu/GeneMark/license_download.cgi"
echo "the gm_key is necessary" 
echo "please download the the GeneMark-ES/ET/EP under the appropriate linux kernel"
echo "it must be copied to your home"
echo -e "do the following:\ngunzip gm_key_64.gz\mv gm_key_64 ~/.gm_key "

#**TSEBRA** available [here](https://github.com/Gaius-Augustus/TSEBRA)
#direct install: 
command='tsebra.py'
if ! command -v $command &> /dev/null
    then
    echo "$command could not be found"
    echo "will try a manual installation through git"
    git clone https://github.com/Gaius-Augustus/TSEBRA 
    cd TSEBRA/bin
    path=$(pwd)
    echo -e "\n#Path to TSEBRA\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../../
fi 

#**trimmomatic** software available [here](http://www.usadellab.org/cms/?page=trimmomatic)
#direct install:
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip


#**GSNAP** for read mappig available [here](http://research-pub.gene.com/gmap/)
#direct install :
command='gsnap'
if ! command -v $command &> /dev/null
    then
    echo "$command could not be found"
    echo "will try a manual installation through git"
    wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2023-10-10.v2.tar.gz
    tar zxf gmap-gsnap-2023-10-10.v2.tar.gz
    cd gmap-2023-10-10/
    path=$(pwd)
    ./configure --prefix=$path  
    make -j8
    make install
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
	cd bin/
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../../
    else
       echo installation failed\nmake sur to have git, make and gclib
       exit 1
    fi
fi

#**gffread** to extract CDS from fasta [see](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
#direct install: 
command='gffread'
if ! command -v $command &> /dev/null
    then
    echo "$command could not be found"
    echo "will try a manual installation through git"
    git clone https://github.com/gpertea/gffread
    cd gffread
    make release
    #if command was successfull then add to path:
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../

    else
       echo installation failed\nmake sur to have git, make and gclib
       exit 1
    fi
fi


## test if orthoFinder is installed
command='orthofinder'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through wget"
    wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder.tar.gz
    tar zxf OrthoFinder.tar.gz
    cd OrthoFinder
    #if command was successfull then add to path:
    path=$(pwd)
    echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    
    #also add diamond, fastme and mcl which are present within orthofinder:
    cd bin/
    path=$(pwd)
    echo -e "\n#Path to diamond\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../../
    #exit 1
fi



#**MCScanX** : 
command='MCScanX'
if ! command -v $command &> /dev/null
then
   #direct install: 
   git clone https://github.com/wyp1125/MCScanX
   cd MCScanX ; make -j 8
   #if command was successfull then add to path:
   if [ $? -eq 0 ]; then
        echo $command installation worked successfully
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../

    else
       echo installation failed\nmake sur to have make and git
       exit 1
    fi
fi


#**muscle** 
command='muscle'
if ! command -v $command &> /dev/null
then
   #direct install: 
   mkdir muscle ; cd muscle
   wget wget https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.linux_intel64
   ln -s muscle5.1.linux_intel64 muscle
   chmod +x muscle
   path=$(pwd)
   echo -e "\n#Path to $command\n export PATH=\$PATH:$path" >> ~/.bashrc 
   source ~/.bashrc  
   cd ../
fi

#**yn00/paml:** 
command='yn00'
if ! command -v $command &> /dev/null
then
   #direct install: 
   wget https://github.com/abacus-gene/paml/releases/download/4.10.7/paml-4.10.7-linux-X86_64.tgz
   tar zxvf paml-4.10.7-linux-X86_64.tgz
   cd paml-4.10.7/bin
   path=$(pwd)
   echo -e "\n#Path to $command\n export PATH=\$PATH:$path" >> ~/.bashrc 
   source ~/.bashrc  
   cd ../../
fi

#**translatorx_vLocal.pl** 
command='translatorx_vLocal.pl'
if ! command -v $command &> /dev/null
then
   #direct install: 
    wget http://161.111.160.230/cgi-bin/translatorx_vLocal.pl
    path=$(pwd)
    echo -e "\n#Path to translatorx\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
fi

#**minimap2** 
#direct install: 
command='minimap2'
if ! command -v $command &> /dev/null
then
   #direct install: 
   curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
   cd ./minimap2-2.26_x64-linux 
   path=$(pwd)
   echo -e "\n#Path to minimap2\n export PATH=\$PATH:$path" >> ~/.bashrc 
   source ~/.bashrc  
   cd ../
fi



#**BUSCO** for quality assesment (https://busco.ezlab.org/)
command='busco'
if ! command -v $command &> /dev/null
then
   #direct install - a new env must be created: 
   mamba create -p busco_env -c bioconda busco=5.5.0
fi 

cd ../

# ----- part for RepeatModeller and RepeatMasker to be completed later ------#
#repeatModeller
git clone https://github.com/Dfam-consortium/RepeatModeler
#install all deps:

wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
tar zxvf RECON-1.08.tar.gz
cd RECON-1.08/src
make 
make install
cd ../bin
if [ $? -eq 0 ]; then
    echo recon installation worked successfully
    path=$(pwd)
    echo -e "\n#Path to Recon\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../
else
   echo installation failed\nmake sur to have make and wget
   exit 1
fi

#repeatscout:
wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
tar zxvf RepeatScout-1.0.6.tar.gz
make
if [ $? -eq 0 ]; then
    echo repeatscout installation worked successfully
    path=$(pwd)
    echo -e "\n#Path to repeatscout\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../
else
   echo installation failed\nmake sur to have make and wget
   exit 1
fi

#trf:
mkdir trf 
cd trf
wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
ln -s trf409.linux64 trf
chmod +x trf
path=$(pwd)
echo -e "\n#Path to trf\nexport PATH=\$PATH:$path" >> ~/.bashrc 
source ~/.bashrc  
#other versions of trf: https://github.com/Benson-Genomics-Lab/TRF/releases/tag/v4.09.1
cd ../

#repeatMasker:
#wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.5.tar.gz

#repeatModeller & repeatMasker dependancies: 
#see : 


#**genemark** 
#wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_pxuuc/gmes_linux_64.tar.gz
echo "to get genemark to work you must register online at http://exon.gatech.edu/GeneMark/license_download.cgi"
echo "the gm_key is necessary" 
echo "please download the the GeneMark-ES/ET/EP under the appropriate linux kernel"
echo "it must be copied to your home"
echo -e "do the following:\ngunzip gm_key_64.gz\mv gm_key_64 ~/.gm_key "

