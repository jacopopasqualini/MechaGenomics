import os
from ftplib import FTP
import subprocess
import argparse

# MAKE DIRS

RES_DIR=['hg37BowtieDB',  'kaijuCore',  'kaijuRefSeq',  'kaijuUHGG',  'krakenRefSeq',  'krakenUHGG',  'taxonkit']

for D in RES_DIR:

    if not os.path.exists(D):

        subprocess.run(['mkdir',D])
        
# CORE KAIJU GETTER

