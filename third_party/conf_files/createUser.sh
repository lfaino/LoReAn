name=$1
uid=$2
adduser --disabled-password --uid $uid --gecos '' $name && adduser $name sudo && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers 
cp /data/gm_key /home/$name/.gm_key
cat /home/$name/.bashrc /etc/profile.d/pathToExport.sh  > /home/$name/.bashrc_new && mv /home/$name/.bashrc_new /home/$name/.bashrc 
su $name
