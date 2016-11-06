#### LOREAN INSTALLATION

### IMPORTANT
LoReAn uses GeneMark-ES as ab-initio software which needs a license key to run. Therefore, IT IS MANDATORY TO download the 64 bit Linux version key for "GeneMark-ES / ET v.4.32" website (http://exon.gatech.edu/GeneMark/license_download.cgi), un-gunzip the key and place it in the folder together with your data.

The best way to use LoReAn is by installing and running the software by Docker
We used Docker because the pipeline uses a lot of software which maybe difficult to install independently.

To install Docker, please refer to:
https://docs.docker.com/engine/installation/

After Docker installation, you can download LoReAn by using:
docker pull lfaino/lorean

Next create a folder in your preferred location (be sure that you have enogh storage space on your drive) and place all your files (short reads, long reads, protein sequence, genome sequences) in the folder. 


Subsequently, from the folder where you input files are, LoReAn can be launched using:
docker run -it --rm -v $PWD:/data lfaino/lorean bash

once inside the container, you can check if LoReAn works using:

lorean.py -help


### KNOWN PROBLEMS 

GMAP compiling
Docker builds container every time the code is updated or modified. Therefore, all the software are compiled using the infrastructure provided by docker. However, the compilation of GMAP is infrastructure dependent. Therefore, if you get a PASA error, try to re-compile GMAP from inside the container using:

cd ~/bin/LoReAn/third_party/software/gmap; make clean; ./configure ; make ; sudo make install


