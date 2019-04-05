# LOREAN INSTALLATION

## IMPORTANT
LoReAn uses GeneMark-ES as *ab-initio* software which needs a license key to run. 

Therefore, **IT IS MANDATORY TO** download the 64 bit Linux version key for [**GeneMark-ES/ET v.4.33 website**](http://exon.gatech.edu/GeneMark/license_download.cgi), un-gunzip the key and place it in the right location.


# LoReAn using SINGULARITY (v2.6.0 - SINGULARITY 3.0-alfa does not work).

### ***We prefer **SINGULARITY** to **DOCKER** because root access is not required.*** 

The best way to use LoReAn is by installing and running the software via [**SINGULARITY**](https://www.sylabs.io/). 
We advice to use **LoReAn** via **SINGULARITY** because the pipeline uses a lot of software which maybe time take to 
install separately. 

A dedicated MYSQL user is required and a linux user is advice. MYSQL user is used by PASA while the Linux user 
is important to not mess-up installations (few files are modified permanently by **SINGULARITY**)

## Here the required steps before using **LoReAn**:



### 2) PLACE THE GENEMARK-ES KEY AT THE RIGHT PLACE 

The next step is to place the ***GeneMark key*** in the home directory of the user running **SINGULARITY**. In Ubuntu, 
the end result would be **~/.gm_key**

If a lorean user is created, add the unzipped gm_key to **/home/lorean** and run the **SINGULARITY** script from the 
lorean home directory (/home/lorean/). 
   

```bash
cd Downloads
gunzip gm_key_64.gz
mv gm_key_64 ~/.gm_key
```

### 3) CREATE MYSQL DATABASE AND DOWNLOAD LOREAN USING SINGULARITY 

These commands can be run from the home directory. The following BASH script will start a new instance of **MYSQL**, download LoReAn
singularity image and move important files.

```bash
singularity pull --bind ${HOME}/mysql/run/mysqld:/run/mysqld/  docker://lfaino/lorean

```



### 4) CHECK THAT LOREAN WORKS

Now, check if  **LoReAn** works
 
 ```bash
lorean.py -help
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

Therefore, **IT IS MANDATORY TO download the 64 bit Linux version key for "GeneMark-ES / ET v.4.33"** website (http://exon.gatech.edu/GeneMark/license_download.cgi), un-gunzip the key and place it in 
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
lorean.py -help
```