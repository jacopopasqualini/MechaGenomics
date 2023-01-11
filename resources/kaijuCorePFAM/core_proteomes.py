import pandas as pd
import ftplib 
from ftplib import FTP
import os
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time

# SERVER
UNIPROT_SITE = 'ftp.uniprot.org'
BAC_DIR = 'pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/'

# LOCAL
UNI_PROTEOMES_DIR='./uniprot_ref_proteomes'

# FIND CORE PFAM
print(">>> Loading PFAM proteomes matrix")
proteomes = pd.read_csv('PFAM_proteome_matrix.csv',index_col='PF')

# GET OCCUPANCY
om = (1*proteomes.fillna(0).astype(bool)).mean(axis=1)
core = list( om[om>0.90].index )
core_proteome = proteomes.loc[core].fillna(0)
roof = core_proteome.max(axis=1)
# SET A MAXIMUM OCCURRENCE, same of the paper
core_families = list(roof[roof<4].index)

# WRITE THE CORE PFAM FAMILIES
print(">>> Core families ")
cpfam_file=open('./sheets/core_pfam.txt','w')
for item in core_families: cpfam_file.writelines(item+'\n')
cpfam_file.close()

i=0
for p in core_families: 
    print(f"  > {i} < {p}")
    i+=1

#anna_families = ['PF00453','PF00572','PF01029','PF01196','PF01649','PF01795','PF03947','PF08338','PF09285','PF17136']

# CROSS INFORMATION BETWEEN PFAM AND UNIPROT GET A LIST OF PROTEOMES
pfam_taxon = list(proteomes.columns)
pfam_taxon = [t.replace('.gz','') for t in pfam_taxon]

uniprot_metadata = pd.read_csv('uniprot_readme_fixed.csv',sep='\t')
uniprot_metadata['Tax_ID']=uniprot_metadata['Tax_ID'].astype(str)
uniprot_metadata=uniprot_metadata.set_index('Tax_ID')

uniprot_metadata = uniprot_metadata.loc[ list(set(pfam_taxon).intersection(uniprot_metadata.index)) ]
uniprot_ids = uniprot_metadata['Proteome_ID'].values

missed_proteomes = set(pfam_taxon).difference(uniprot_metadata.index)

#missed_pfam_proteomes_file=open('./sheets/missed_pfam_proteomes.txt','w')
#for item in missed_proteomes:  missed_pfam_proteomes_file.writelines(item+'\n')
#missed_pfam_proteomes_file.close()

with open('./sheets/missed_pfam_proteomes.txt','w') as mpp_f:
    for item in missed_proteomes:  mpp_f.writelines(item+'\n')
mpp_f.close()

##### SCRAPING UNIPROT

ftp = FTP(UNIPROT_SITE)
ftp.login()

pid = uniprot_metadata['Proteome_ID']

core_uniproteomes={}

multifasta_file=open('./sheets/multiple_reference.txt','w')
corepfams_reduced=open('corepfams_reduced.fa','w')
empty_proteomes=open('./sheets/empty_proteomes.txt','w')
redundant_prots=open('./sheets/redundant.txt','w')


for b in range(len(uniprot_ids)):

    # track the process
    p = uniprot_ids[b]
    t = list(pid[pid==p].index)[0]
    print(f"  > {b} < uniprotId: {p} < taxonId: {t}")
    
    ftp.cwd(BAC_DIR+p)

    # find proteome.fasta
    files = []
    try: files = ftp.nlst()
    except ftplib.error_perm as resp:
        if str(resp) == "550 No files found": print("No files in this directory")
        else: raise

    fasta = [f for f in files if (('fasta' in f) and ('DNA'  not in f))]

    # report if multiple fasta are available
    if len(fasta)>1: 
        ff = ''
        for f in fasta: ff+=(f+'\t')
        multifasta_file.writelines(f'{p}\t{ff}\n')

    proteome_fasta = fasta[0]
    
    # download proteome.fasta
    if  not os.path.exists( proteome_fasta ):

        while not os.path.exists( proteome_fasta ):

            try: os.system(f'wget https://{UNIPROT_SITE}/{BAC_DIR}/{p}/{proteome_fasta}')
            except TimeoutError: time.sleep(1)

    os.system(f'mv {proteome_fasta} {UNI_PROTEOMES_DIR}/{proteome_fasta} ')
    ftp.cwd('/')

    # pick up the corresponding proteome with same ncbi numerical taxon id
    # keep only its core-pfam part
    # such PFAMs will appear in specific proteins, we will look for such proteins in the uniprot proteome
    # in the prfam proteome the sequence position is reported: the cp (core pfam) adress will save such position for each pfam in each protein
 
    pfam_proteome = pd.read_csv(f'./proteomes/{t}.tsv.gz',compression='gzip',sep='\t',skiprows=3,header=None)
    
    core_proteome = pfam_proteome[pfam_proteome[5].isin(core_families)]

    #core_pfam = core_proteome[[0,5]]
    #core_pfam = core_pfam.set_index(0)
    #core_pfam = core_pfam.groupby(level=0).agg(list).to_dict()
    
    proteins = core_proteome[0].unique()
    
    cp_adress = {} 

    for pc in proteins:

        A=core_proteome[ core_proteome[0]==pc ]
        adress = {}

        for a in A.index:
            ad = []
            ad.append(A.loc[a,1])
            ad.append(A.loc[a,2])
            pf = A.loc[a,5]
            adress[pf]=ad

        cp_adress[pc]=adress

    #print(cp_adress)
    # open proteome fasta file with biopython. 
    # since bacterial proteomes are not that big, we can store them in vocabularies
    # for each taxon t save in the core_uniproteomes the corresponding core sequences
    
    
    with gzip.open(f'{UNI_PROTEOMES_DIR}/{proteome_fasta}', "rt") as fastagz:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fastagz, "fasta"))
        fasta_ids = list(fasta_dict.keys())


    for p in cp_adress.keys():

        in_fasta_proteins = [f for f in fasta_ids if p in f]

        if len(in_fasta_proteins)>1: redundant_prots.write(f'{t} \t {in_fasta_proteins}')
    
        for f in cp_adress[p].keys():            

            a = cp_adress[p][f]
            try:

                core_seq = fasta_dict[in_fasta_proteins[0]].seq[a[0]:a[1]]
                rec = SeqRecord(seq=core_seq, id=f'{p}.{f}_{t}', name='', description='', dbxrefs=[])
                #print(rec)
                SeqIO.write(rec,corepfams_reduced,"fasta")

            except IndexError:
                print('Missed protein/proteome')
                empty_proteomes.writelines(item+'\n')
            

corepfams_reduced.close()
empty_proteomes.close()
redundant_prots.close()