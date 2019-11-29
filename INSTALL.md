# LOREAN INSTALLATION

## IMPORTANT
LoReAn uses GeneMark-ES as *ab-initio* software which needs a license key to run. 

Therefore, **IT IS MANDATORY TO** download the 64 bit Linux version key for [**GeneMark-ES/ET v.4.33 website**](http://exon.gatech.edu/GeneMark/license_download.cgi), un-gunzip the key and place it in the right location.


# LoReAn using SINGULARITY.

### ***We prefer **SINGULARITY** to **DOCKER** because root access is not required.*** 

The best way to use LoReAn is by installing and running the software via [**SINGULARITY**](https://www.sylabs.io/). 
We advice to use **LoReAn** via **SINGULARITY** because the pipeline uses a lot of software which maybe time take to 
install separately. 

## Here the required steps before using **LoReAn**:

### 1) PLACE THE GENEMARK-ES KEY AT THE RIGHT PLACE 

The first step is to place the ***GeneMark key*** in the home directory of the user running **SINGULARITY**. In Ubuntu, 
the end result would be **~/.gm_key**
   

```bash
cd Downloads
gunzip gm_key_64.gz
mv gm_key_64 ~/.gm_key
```

### 2) DOWNLOAD THE REQUIRED FILES/FOLDERS

In order to run **LoReAn** by **Singularity exec** command, you need to download and unzip these two files:

```bash
wget https://github.com/lfaino/LoReAn/raw/noIPRS/third_party/software/config.augustus.tar.gz && tar -zxvf config.augustus.tar.gz
```
```bash
wget https://github.com/lfaino/LoReAn/raw/noIPRS/third_party/software/RepeatMasker.Libraries.tar.gz && tar -zxvf RepeatMasker.Libraries.tar.gz
```    
The firts file is the configuration folder from Augustus software (see below <PATH_TO_AUGUSTUS_CONF_FOLDER>) while the 
second file is the Libraries folder of RepeatMasker software (see below <PATH_TO_LIBRARY_FOLDER>)

Next you can download and build the **Syngularity** image using:  

```bash
singularity pull docker://lfaino/lorean:latest
```


### 3) CHECK THAT LOREAN WORKS

Now, check if  **LoReAn** works by
 
 ```bash
 
singularity exec -B <PATH_TO_AUGUSTUS_CONF_FOLDER>:/opt/LoReAn/third_party/software/augustus/config/ -B 
<PATH_TO_LIBRARY_FOLDER>:/usr/local/RepeatMasker/Libraries/ <PATH_TO_LOREAN_IMAGE>/lorean_latest.sif lorean -h

 ```

At this point, you should see the options list. 
You can continue by testing lorean using the toy datasets located at [LoReAn examples](https://github.com/lfaino/LoReAn_Example)



## LoReAn using Docker.

On Linux system, make sure that the user runnig docker is added as user in docker user group.

On Windows system, before installing Docker **IT IS MANDATORY** to allow symbolic links. PASA makes symbolic during the run.
The eisiest way to run docker is via Docker Toolbox. During the installation, set the size of the disk image to about 30Gb.
After the installation run Docker Quickstart Terminal and follow the instruction below 

### IMPORTANT
LoReAn uses GeneMark-ES as ab-initio software which needs a license key to run. 

Therefore, **IT IS MANDATORY TO download the 64 bit Linux version key for "GeneMark-ES / ET v.4.48"** website (http://exon.gatech.edu/GeneMark/license_download.cgi), un-gunzip the key and place it in 
the folder together with your data.

The best way to use LoReAn is by installing and running the software by Docker
We used Docker because the pipeline uses a lot of software which maybe difficult to install independently.


To install Docker, please refer to:
https://docs.docker.com/engine/installation/

After Docker installation, you can download  LoReAn by using:
```bash
docker run -it --rm -v $PWD:/data lfaino/lorean createUser.py $USER $UID 
```

At this point, run

```bash
lorean -help
```