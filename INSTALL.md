# LOREAN INSTALLATION

## IMPORTANT
LoReAn uses GeneMark-ES as *ab-initio* software which needs a license key to run. 

Therefore, **IT IS MANDATORY TO download the 64 bit Linux version key for "GeneMark-ES / ET v.4.33"** website 
(http://exon.gatech.edu/GeneMark/license_download.cgi), un-gunzip the key and place it in the right location.


# LoReAn using Singularity.

The best way to use LoReAn is by installing and running the software by **SINGULARITY**. 
We advice to use **LoReAn** via **SINGULARITY** because the pipeline uses a lot of software which maybe difficult to 
install all of them independently. We prefer **SINGULARITY** to **DOCKER** because root access is not required.

However, a dedicated MYSQL user is required and a linux user is advice. Few steps are required before using **LoReAn**
MYSQL user is used by PASA while the Linux user is important to not mess-up installations

MYSQL user:
```bash
CREATE USER 'lorean'@'localhost' IDENTIFIED BY 'lorean';
GRANT ALL PRIVILEGES ON * . * TO 'lorean'@'localhost';
FLUSH PRIVILEGES;
```
   
After **MYSLQ** user is created, we can start **SYNGULARITY**.

```bash
singularity shell --bind /var/run/mysqld/:/run/mysqld/  docker://lfaino/lorean:iprscan_rpMask
```

At this point, some files need to be moved
```bash
cat /home/lorean/.bashrc /etc/profile.d/pathToExport.sh  > /home/lorean/.bashrc_new && mv /home/lorean/.bashrc_new /home/lorean/.bashrc
source ~/.bashrc
cp -r /opt/LoReAn/third_party/software/augustus/ /home/lorean/
```
 




# LoReAn using Docker.

### On Linux

On Linux system, make sure that the user is added as user in docker user group.

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
if panther is not present. **Panther** is a big database and the user need to download it 

At this point, run

```bash
lorean.py -help
```

###IMPORTANT
Before running the docker command, place the unzipped GeneMark key in the folder from where you are running the docker 
command 

###On Windows
On Windows system, before installing Docker **IT IS MANDATORY** to allow symbolic links. PASA makes symbolic during the run.
The esiest way to run docker is via Docker Toolbox. During the installation, set the size of the disk image to about 30Gb.
After the installation run Docker Quickstart Terminal and follow the instruction to run on Linux (above)

###IMPORTANT
Before running the docker command, place the unzipped GeneMark key in the folder from where you are running the docker 
command 


## Run LoReAn in standard bash

It is possible to run LoReAn using bash only if all the tools listed in the Dockerfile are installed, the specific file removed and
each tool is correctly installed. After the installation is complete, LoReAn can be run  


