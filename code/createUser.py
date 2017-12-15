#! /usr/bin/env python3


import os
import shutil
import subprocess
import sys
from pathlib import Path


###############
###FUNCTIONS###
###############


def create_user():

    name_user = sys.argv[1]
    uid_user = sys.argv[2]

    root = os.getcwd()
    sys.stdout.write(('\n### CREATING USER WITH NAME %s AND UID %s IN THE DOCKER IMAGE ###\n\n') % (name_user, uid_user))

    log_file = os.path.join(root, "CreateUser.log.txt")
    err_file = os.path.join(root, "CreateUser.err.txt")
    log = open(log_file, 'w')
    err = open(err_file, 'w')

    com = "adduser --disabled-password --uid %s --gecos '' %s && adduser %s sudo && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers" % (uid_user, name_user, name_user, "%s")
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "cp /data/gm_key /home/%s/.gm_key" % (name_user)
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    gm_key_file = Path("/home/%s/.gm_key" % (name_user))
    if not gm_key_file.is_file():
        sys.exit("#####PLEASE PLACE THE gm_key IN THE DIRECTORY WITH ALL THE OTHER FILES.#####\n")

    com = "cat /home/%s/.bashrc /etc/profile.d/pathToExport.sh  > /home/%s/.bashrc_new && mv /home/%s/.bashrc_new /home/%s/.bashrc" % (name_user, name_user, name_user, name_user)
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "usermod -d /var/lib/mysql/ mysql"
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "/etc/init.d/mysql start"
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "mysql --user=\"root\" --password=\"lorean\" --execute=\"set global sql_mode='STRICT_TRANS_TABLES,ERROR_FOR_DIVISION_BY_ZERO,NO_AUTO_CREATE_USER,NO_ENGINE_SUBSTITUTION';\""
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "chown -R %s:%s /home/%s/.gm_key" % (name_user, name_user, name_user)
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "chown -R %s:%s /home/%s/.bashrc" % (name_user, name_user, name_user)
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "cp -r /opt/LoReAn/third_party/software/augustus/ /home/%s/" % (name_user)
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "chmod -R 775 /opt/LoReAn"
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "chmod -R 775 /home/%s/augustus" % (name_user)
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    com = "chown -R %s:%s /home/%s/augustus" % (name_user, name_user, name_user)
    create_user_call = subprocess.Popen(com, stdout=log, stderr=err, shell=True)
    create_user_call.communicate()

    shutil.chown(log_file, user=name_user, group=name_user)
    shutil.chown(err_file, user=name_user, group=name_user)

    subprocess.run(["su", name_user])


if __name__ == '__main__':
    create_user()