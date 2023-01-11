from socketserver import ThreadingUnixStreamServer
from turtle import back
import pytaxonkit

import pandas as pd
import os
import time
import ftplib 
from ftplib import FTP

import subprocess


LOCAL_DIR = './'
FTP_SITE  = 'ftp.ebi.ac.uk'
FTP_PROTEOMES = '/pub/databases/Pfam/current_release/proteomes'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# connect with ftp pfam database and go to the proteomes folder

print(">>> PFAM-FTP SITE LOG-IN: ")
ftp = FTP(FTP_SITE)
ftp.login()
ftp.cwd(FTP_PROTEOMES)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# get all the taxa in the protomes folder
print(">>> PICKING BACTERIAL PROTEOMES")

files_list = ftp.nlst()

#print(files_list)

taxa_list = [ f[:-len('.tsv.gz')] for f in files_list ]

del files_list
# get the ncbi taxonomy table
print(">>> NCBI TAXONOMY FILTERING")

TAXON_DIR = '/mnt/MetaGym/ecogut/resources/taxonkit'
taxa_atlas=pytaxonkit.lineage(ids=taxa_list,data_dir=TAXON_DIR)
# keep only the taxonomic path and split it 
taxa_atlas=taxa_atlas.set_index('TaxID')

lineage = taxa_atlas['Lineage'].str.contains('Bacteria')
taxonk_bacteria = list( lineage[lineage==True].index )


#DIR = '/Users/jpasqual/Downloads/'
print(">>> UNIPROT REFERENCE PROTEOMES")

uniprot_file = os.path.join( LOCAL_DIR,'uniprot_reference_proteomes.csv')
uniprot=pd.read_csv(uniprot_file,sep='\t')
nr_uniprot_bacteria=set(uniprot['Taxon'].values.astype('str'))

print()
print('  > PFAM proteomes:',len(taxa_list))
print('  > nr uniprot proteomes:',len(nr_uniprot_bacteria))
print('  > Bacterial PFAM proteomes with taxonkit: ',len(taxonk_bacteria))
print('  > Non redundant proteomes in PFAM',len(list(set(nr_uniprot_bacteria).intersection(taxa_list))))

pfam_bacteria = list(set(nr_uniprot_bacteria).intersection(taxonk_bacteria))
print('  > Non redundant proteomes in PFAM with taxonkit',len(pfam_bacteria))
print()

#write all bacterial taxa into a file
#textfile = open(LOCAL_DIR+'/pfam_bacteria.txt', 'w')
#for b in pfam_bacteria:
#    textfile.write(b + "\n")
#textfile.close()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# download all bacterial proteomes

PROTEOMES_DIR  = LOCAL_DIR +'proteomes/'
PROTEOMES_PAGE = 'http://' + FTP_SITE + FTP_PROTEOMES + '/'
Abundances = pd.DataFrame()

missed = []

B=len(pfam_bacteria)

for b in range(B):
    
    p=pfam_bacteria[b]
    
    proteome_ftp_file = PROTEOMES_PAGE + p + '.tsv.gz'
    proteome_loc_file = PROTEOMES_DIR  + p + '.tsv.gz'
    
    print('\n','* '*30)
    print(proteome_loc_file)
    print(b,')',p)
    
    # download the file
    while not os.path.exists( proteome_loc_file ):

        os.system(f'wget {proteome_ftp_file}')
        os.system(f'mv {p}.tsv.gz ./proteomes/{p}.tsv.gz')