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


To install Docker, please refer to:
https://docs.docker.com/engine/installation/

After Docker installation, you can download  LoReAn by using:
```bash
docker run -it --rm -v $PWD:/data -v /path/to/panther/folder/panther:/data_panther lfaino/lorean:iprscan_rpMask createUser.py $USER $UID
```

or
```bash
docker run -it --rm -v $PWD:/data lfaino/lorean:iprscan_rpMask createUser.py $USER $UID
```
if panther is not present

At this point, run

```bash
lorean.py -help
```

## Run LoReAn in standard bash

It is possible to run LoReAn using bash only if all the tools listed in the Dockerfile are installed, the specific file removed and
each tool is correctly installed. After the installation is complete, LoReAn can be run  

