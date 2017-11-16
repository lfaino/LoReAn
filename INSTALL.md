# LOREAN INSTALLATION

## Use LoReAn using Docker.

On Linux system, make sure that the user runnig docker is added as user in docker user group.

On Windows system, before installing Docker **IT IS MANDATORY** to allow symbolic links. PASA makes symbolic during the run.
The eisiest way to run docker is via Docker Toolbox. During the installation, set the size of the disk image to about 30Gb.
After the installation run Docker Quickstart Terminal and follow the instruction below 

### IMPORTANT
LoReAn uses GeneMark-ES as ab-initio software which needs a license key to run. 

Therefore, **IT IS MANDATORY TO download the 64 bit Linux version key for "GeneMark-ES / ET v.4.33"** website (http://exon.gatech.edu/GeneMark/license_download.cgi), un-gunzip the key and place it in 
the folder together with your data.

The best way to use LoReAn is by installing and running the software by Docker
We used Docker because the pipeline uses a lot of software which maybe difficult to install independently.

Before installing Docker, please check your UID:
```bash
id user_name
```
where user_name is the name of the user that runs LoReAn on the host machine

To install Docker, please refer to:
https://docs.docker.com/engine/installation/

After Docker installation, you can download and run LoReAn by using:
```bash
docker run -it --rm -v $PWD:/data lfaino/lorean createUser.sh user_name uid_user
```

At this point, run

```bash
lorean.py -help
```

## Run LoReAn in standard bash




### KNOWN PROBLEMS 

GMAP compiling:

Docker builds container every time the code is updated or modified. Therefore, all the software are compiled using the infrastructure provided by docker. However, the compilation of GMAP is 
infrastructure dependent. Therefore, if you get a PASA error and the file **annotation/run/PASA/gmap.spliced_alignments.gff3** is empty, try to re-compile GMAP from inside the 
container using:
```bash
exit
cd /opt/LoReAn/third_party/software/gmap; make clean; ./configure ; make ; sudo make install ; cd /data/
su user_name
```

