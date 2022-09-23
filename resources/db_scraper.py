import os
from ftplib import FTP
import subprocess
import argparse

# MAKE DIRS

RES_DIR=['hg37BowtieDB',  'kaijuCore',  'kaijuRefSeq',  'kaijuUHGG',  'krakenRefSeq',  'krakenUHGG',  'taxonkit']

for D in RES_DIR:

    if not os.path.exists(D):

        subprocess.run(['mkdir',D])

# REFSEQ KAIJU GETTER
# https://kaiju.binf.ku.dk/server
         
kaijuRefSeq_db_page = 'https://kaiju.binf.ku.dk/database/kaiju_db_refseq_2022-03-23.tgz'

kaijuRefSeq_files=['kaiju_db_refseq.fmi', 'names.dmp', 'nodes.dmp']

if os.listdir('kaijuRefSeq')!=kaijuRefSeq_files:

    subprocess.run(['wget','https://kaiju.binf.ku.dk/database/kaiju_db_refseq_2022-03-23.tgz'])

# CORE KAIJU GETTER



