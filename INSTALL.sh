### Install openjdk.
sudo apt install openjdk-8-jdk # can possibly work with higher versions

### Python modules.
sudo apt install python3-pip python3-dev
sudo -H python3 -m pip install pysam pandas biopython

### Install third-party software
sudo apt install bwa samtools  

### Unzip reference
sourcedir=$(dirname $0)
gunzip ${sourcedir}/reference/*
