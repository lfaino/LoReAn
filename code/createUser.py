#! /usr/bin/env python3


import os
import subprocess
import sys
import tempfile
from pathlib import Path


###############
###FUNCTIONS###
###############


def create_user():

    name_user = sys.argv[1]
    uid_user = sys.argv[2]

    root = os.getcwd()
    sys.stdout.write(('\n### CREATING USER WITH NAME %s AND UID %s IN THE DOCKER IMAGE ###\n') % (name_user, uid_user))

    if len(sys.argv) > 3 and "--verbose" in sys.argv[3] or len(sys.argv) > 3 and "-v" in sys.argv[3]:
        log = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=root, prefix="startUser.", suffix=".out")
        err = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=root, prefix="startUser.", suffix=".out")
    else:
        log = tempfile.NamedTemporaryFile(mode='w', dir=root, prefix="startUser.", suffix=".out")
        err = tempfile.NamedTemporaryFile(mode='w', dir=root, prefix="startUser.", suffix=".out")

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

    os.system("su %s >> %s" % (name_user, log.name))


if __name__ == '__main__':
    create_user()