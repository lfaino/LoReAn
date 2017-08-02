#!/usr/bin/env bash
name=$1
uid=$2
adduser --disabled-password --uid $uid --gecos '' $name && adduser $name sudo && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers 
cp /data/gm_key /home/$name/.gm_key
cat /home/$name/.bashrc /etc/profile.d/pathToExport.sh  > /home/$name/.bashrc_new && mv /home/$name/.bashrc_new /home/$name/.bashrc 
usermod -d /var/lib/mysql/ mysql
/etc/init.d/mysql start
mysql --user="root" --password="lorean" --execute="set global sql_mode='STRICT_TRANS_TABLES,ERROR_FOR_DIVISION_BY_ZERO,NO_AUTO_CREATE_USER,NO_ENGINE_SUBSTITUTION';"
chown -R $name:$name /home/$name/.gm_key
chown -R $name:$name /home/$name/.bashrc
#chown -R $name:$name /opt/LoReAn/
chmod -R 775 /opt/LoReAn
su $name
