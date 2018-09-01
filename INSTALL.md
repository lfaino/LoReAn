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

The next step is to place the ***GeneMark key*** in the home directory of the user running **SINGULARITY**. If a user lorean is created,
add the unzipped gm_key to **/home/lorean** and run the **SINGULARITY** script from the lorean home directory (/home/lorean/). 
In Ubuntu, the end result would be **/home/lorean/.gm_key**   

```bash
cd Downloads
gunzip gm_key_64.gz
mv gm_key_64 ~/.gm_key
```

### 2) CREATE MYSQL DATABASE AND DOWNLOAD LOREAN USING SINGULARITY 


```bash
singularity pull --name mysql.simg shub://ISU-HPC/mysql
wget -O ./.my.cnf https://raw.githubusercontent.com/lfaino/LoReAn/dev/third_party/conf_files/my.cnf 
wget -O ./.mysqlrootpw https://raw.githubusercontent.com/lfaino/LoReAn/dev/third_party/conf_files/mysqlrootpw
mkdir -p ${PWD}/mysql/var/lib/mysql ${PWD}/mysql/run/mysqld
singularity instance.start --bind ${HOME} --bind ${PWD}/mysql/var/lib/mysql/:/var/lib/mysql --bind ${PWD}/mysql/run/mysqld:/run/mysqld ./mysql.simg mysql
singularity run instance://mysql
singularity shell --bind ${PWD}/mysql/run/mysqld:/run/mysqld/  docker://lfaino/lorean:iprscan_rpMask
cat /home/lorean/.bashrc /opt/LoReAn/third_party/conf_files/pathToExport.txt  > /home/lorean/.bashrc.lorean
source ~/.bashrc.lorean
cp -r /opt/LoReAn/third_party/software/augustus/ /home/lorean/
```

After the first use:
```bash
source ~/.bashrc.lorean
```

### 3) CHECK THAT LOREAN WORKS

Now, check if  **LoReAn** works
 
 ```bash
lorean.py -help
 ```

At this point, you should see the options list. 
You can continue by testing lorean using the toy datasets located at [LoReAn examples](https://github.com/lfaino/LoReAn_Example)
