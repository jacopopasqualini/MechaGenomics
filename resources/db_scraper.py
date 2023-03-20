import os
from ftplib import FTP
import subprocess
import argparse

# MAKE DIRS

RES_DIR=['hg37BowtieDB',  'kaijuCore',  'kaijuRefSeq',  'kaijuUHGG',  'krakenRefSeq',  'krakenUHGG']

for D in RES_DIR:

    if not os.path.exists(D):

        subprocess.run(['mkdir',D])

# REFSEQ KAIJU GETTER
# info @: https://kaiju.binf.ku.dk/server
        
kaijuRefSeq_files=['kaiju_db_refseq.fmi', 'names.dmp', 'nodes.dmp']

if os.listdir('kaijuRefSeq')!=kaijuRefSeq_files:

    #subprocess.run(['wget','-r','https://kaiju.binf.ku.dk/database/kaiju_db_refseq_2022-03-23.tgz'])

    if os.path.exists('./kaiju_db_refseq_2022-03-23.tgz')
        subprocess.run(['mv','kaiju_db_refseq_2022-03-23.tgz','kaijuRefSeq/kaiju_db_refseq_2022-03-23.tgz'])
        subprocess.run(['tar','-zxvf','kaijuRefSeq/kaiju_db_refseq_2022-03-23.tgz'])

'''
# CORE KAIJU GETTER

#kaijuRefSeq_files=['kaiju_db_refseq.fmi', 'names.dmp', 'nodes.dmp']

if os.listdir('kaijuRefSeq')!=kaijuRefSeq_files:

    subprocess.run(['wget','https://github.com/liphlab/Kaiju-core/archive/refs/heads/master.zip'])

# KRAKEN2 STANDARD

condition=False

if condition==True:

    subprocess.run(['wget','https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220908.tar.gz'])


# KRAKEN UHGG

if condition==True:

    subprocess.run(['wget','http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/kraken2_db_uhgg_v2/']=

'''