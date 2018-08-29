# LOREAN INSTALLATION

## IMPORTANT
LoReAn uses GeneMark-ES as *ab-initio* software which needs a license key to run. 

Therefore, **IT IS MANDATORY TO download the 64 bit Linux version key for "GeneMark-ES / ET v.4.33"** website 
(http://exon.gatech.edu/GeneMark/license_download.cgi), un-gunzip the key and place it in the right location.


# LoReAn using Singularity (v2.6.0).

### ***We prefer **SINGULARITY** to **DOCKER** because root access is not required.*** 

The best way to use LoReAn is by installing and running the software via **SINGULARITY** (https://www.sylabs.io/). 
We advice to use **LoReAn** via **SINGULARITY** because the pipeline uses a lot of software which maybe time take to 
install separately. 

A dedicated MYSQL user is required and a linux user is advice. MYSQL user is used by PASA while the Linux user 
is important to not mess-up installations (few files are modified permanently by **SINGULARITY**)

## Here the required steps before using **LoReAn**:

### 1a) CREATE MYSQL DATABASE 

How to create a **MYSQL** user in the host system:
```bash
CREATE USER 'lorean'@'localhost' IDENTIFIED BY 'lorean';
GRANT ALL PRIVILEGES ON * . * TO 'lorean'@'localhost';
FLUSH PRIVILEGES;
```

### 1b) CREATE LOREAN USER (UBUNTU)(OPTIONAL - not sure if possible in "BASH on Ubuntu on Windows 10") 

```bash
sudo adduser lorean
```

After **MYSLQ** and **Linux** user are created, we can start **SYNGULARITY**. 

### 2) PLACE THE GENEMARK-ES KEY AT THE RIGHT PLACE 

The next step is to place the ***GeneMark key*** in the home directory of the user running **SINGULARITY**. If a user lorean is created,
add the unzipped gm_key to **/home/lorean** and run the **SINGULARITY** script from the lorean home directory (/home/lorean/). 
In Ubuntu, the end result would be **/home/lorean/.gm_key**   

```bash
cd Downloads
gunzip gm_key_64.gz
mv gm_key_64 ~/.gm_key
```


### 3) DOWNLOAD AND START LOREAN SHELL VIA SINGULARITY  

***IT IS MANDATORY*** to bind MYSQL running on the host to the **SYNGULARITY** image. To do so, search for the **mysqld.sock** file
(In UBUNTU is located at /run/mysqld/). Use the --bind option to link the folder containing the **mysqld.sock** to the 
image (see command below)

```bash
singularity shell --bind /var/run/mysqld/:/run/mysqld/  docker://lfaino/lorean:iprscan_rpMask
```

### 4) MOVE IMPORTANT FILES 

This step need to be performed only the first time **SINGULARITY** is run. The above changes are stored permanently 
in the lorean home directory and used in all following runs. Therefore, we suggest to have a dedicated home directory 
to run **SINGULARITY** 

At this point, some files need to be moved
```bash
cat /home/lorean/.bashrc /etc/profile.d/pathToExport.sh  > /home/lorean/.bashrc_new
mv /home/lorean/.bashrc /home/lorean/.bashrc.bk
mv /home/lorean/.bashrc_new /home/lorean/.bashrc
source ~/.bashrc
cp -r /opt/LoReAn/third_party/software/augustus/ /home/lorean/
```
  
### 4) CHECK THAT LOREAN WORKS

Now, check if  **LoReAn** works
 
 ```bash
lorean.py -help
 ```

At this point, you should see the options list. 
You can continue by testing lorean using the toy datasets located at: [LoReAn examples](https://github.com/lfaino/LoReAn_Example)
