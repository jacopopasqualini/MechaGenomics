from re import sub
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

LOCAL_PRT_DIR = './species_proteomes'

# fallo con os.sys ecc in modo da automatizzare creazione indice e unisci questa cartella con quella di kaiju.
# sempre in questa carttella scarica la tassonomia di kraken
# scarica in questa cartella il file prelim map ecrea con seqkit la mappa tra codice genoma e tassonomia
# apri un nuovo fastq  dove salvi le proteine con loro id di uhgg e attaccato quello di kraken fatto da prelimminary map

UHGG_FTP="http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/kraken2_db_uhgg_v2/taxonomy/"

if not os.path.exists("prelim_map.csv"):

    MAP_TXT = f"{UHGG_FTP}prelim_map.txt"
    os.system(f"wget {MAP_TXT}")
    os.system("mv prelim_map.txt prelim_map.csv")

if not os.path.exists("nodes.dmp"):

    NODES= f"{UHGG_FTP}nodes.dmp"
    os.system(f"wget {NODES}")

if not os.path.exists("names.dmp"):

    NAMES= f"{UHGG_FTP}names.dmp"
    os.system(f"wget {NAMES}")

kraken_map = pd.read_csv("prelim_map.csv",sep='\t',header=None)

kraken_map[3]=kraken_map[1].str.split("|",expand=True)[0].str.split("_",expand=True)[0]
kraken_map=kraken_map[[2,3]].drop_duplicates()
kraken_map=kraken_map.rename(columns={2:"tax_id",3:"uhgg_id"})
kraken_map=kraken_map.set_index("uhgg_id")
kraken_map=kraken_map['tax_id'].astype(str)
kraken_map=kraken_map.to_dict()

with open('UHGGv2.faa', "w") as uhgg_fasta:

    for f in os.listdir("./species_proteomes"):

        print(f.replace(".faa",""))

        with open(f'./species_proteomes/{f}', "r") as fasta:

            fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
            fasta_ids = list(fasta_dict.keys())

            taxa_id = kraken_map[ f.replace(".faa","") ]
            
            for id in fasta_dict.keys(): 
 
                rec = SeqRecord(seq=fasta_dict[id].seq, id=f'{id.replace("_",".")}_{taxa_id}', name='', description='', dbxrefs=[])
                SeqIO.write(rec,uhgg_fasta,"fasta")

uhgg_fasta.close()