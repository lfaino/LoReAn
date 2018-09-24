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


### 1) CREATE LOREAN USER (UBUNTU) - (THIS IS AN OPTIONAL STEP) 

```bash
sudo adduser lorean
```

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

**NOTE**: **MYSQL** will run on port 5123 to avoid conflict with other **MYSQL** instance already running on the system. Please, 
check that the door is open and available to use

**NOTE**: the original .barsrc file should be without any **export $PATH:** add to it. If you added personal PATH to 
the ~/.bahsrc using the **export** command, please remove them from the final **~/.bashrc.lorean** before running the below 
source command.   

```bash
singularity pull --name mysql.simg shub://ISU-HPC/mysql
wget -O ~/.my.cnf https://raw.githubusercontent.com/lfaino/LoReAn/dev/third_party/conf_files/my.cnf 
wget -O ~/.mysqlrootpw https://raw.githubusercontent.com/lfaino/LoReAn/dev/third_party/conf_files/mysqlrootpw
mkdir -p ${HOME}/mysql/var/lib/mysql ${HOME}/mysql/run/mysqld
singularity instance.start --bind ${HOME} --bind ${HOME}/mysql/var/lib/mysql/:/var/lib/mysql --bind ${HOME}/mysql/run/mysqld:/run/mysqld ./mysql.simg mysql
singularity run instance://mysql
singularity shell --bind ${HOME}/mysql/run/mysqld:/run/mysqld/  docker://lfaino/lorean:iprscan_rpMask
cat ~/.bashrc /opt/LoReAn/third_party/conf_files/pathToExport.txt  > ~/.bashrc.lorean
source ~/.bashrc.lorean
cp -r /opt/LoReAn/third_party/software/augustus/ ~/
```



### 4) CHECK THAT LOREAN WORKS

Now, check if  **LoReAn** works
 
 ```bash
lorean.py -help
 ```

At this point, you should see the options list. 
You can continue by testing lorean using the toy datasets located at [LoReAn examples](https://github.com/lfaino/LoReAn_Example)


### OTHER USEFUL COMMANDS

###TO STOP MYSQL INSTANCE 
Use the following command to stop mysql instance.

```bash
singularity instance.stop mysql
```
### TO START LOREAN AFTER THE FIRST USE

After the first use, the augustus folder and the bashrc.lorean are already prepared. therefore we need only to download 
the singularity images. Here the code:

```bash
singularity pull --name mysql.simg shub://ISU-HPC/mysql
wget -O ~/.my.cnf https://raw.githubusercontent.com/lfaino/LoReAn/dev/third_party/conf_files/my.cnf 
wget -O ~/.mysqlrootpw https://raw.githubusercontent.com/lfaino/LoReAn/dev/third_party/conf_files/mysqlrootpw
mkdir -p ${HOME}/mysql/var/lib/mysql ${HOME}/mysql/run/mysqld
singularity instance.start --bind ${HOME} --bind ${HOME}/mysql/var/lib/mysql/:/var/lib/mysql --bind ${HOME}/mysql/run/mysqld:/run/mysqld ./mysql.simg mysql
singularity run instance://mysql
singularity shell --bind ${HOME}/mysql/run/mysqld:/run/mysqld/  docker://lfaino/lorean:iprscan_rpMask
source ~/.bashrc.lorean
```


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
docker run -it --rm -v $PWD:/data lfaino/lorean:iprscan_rpMask createUser.py $USER $UID 
```

At this point, run

```bash
lorean.py -help
```