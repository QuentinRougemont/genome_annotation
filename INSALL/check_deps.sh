#check that everything is installed!

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

#**samtools**
command='samtools'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi 

#**R** with several packages 
command='R'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi 

#**transeq**  from EMBOSS to convert fasta into protein [click_here](https://www.bioinformatics.nl/cgi-bin/emboss/help/transeq) 
command='transeq'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi
 
command='braker.pl'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi

#**Protint** 
command='prothint.py'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi

#**Diamond**
command='diamond'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi


#**cdbfasta**:
command='cdbfasta'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi

#**bamtools**
command='bamtools'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi


#**Augustus**
command='augustus'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi

#**genemark** 
command="gms2hints.pl" 
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi

#**TSEBRA** available [here](https://github.com/Gaius-Augustus/TSEBRA)
command='tsebra.py'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi 

#**GSNAP** for read mappig available [here](http://research-pub.gene.com/gmap/)
command='gsnap'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi

#**gffread** to extract CDS from fasta [see](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
command='gffread'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi


### test if orthoFinder is installed
command='orthofinder'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi
#
#
#
##**MCScanX** : 
command='MCScanX'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi
#
##**muscle** 
command='muscle'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi
#
##**yn00/paml:** 
command='yn00'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi
#
##**translatorx_vLocal.pl** 
command='translatorx_vLocal.pl'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi
#
##**minimap2** 
##direct install: 
command='minimap2'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi
#
##**BUSCO** for quality assesment (https://busco.ezlab.org/)
command='busco'
if ! command -v $command &> /dev/null
then
    echo "ERROR: $command could not be found"
    echo "please install it before doing anything else"
    exit 1
fi 
