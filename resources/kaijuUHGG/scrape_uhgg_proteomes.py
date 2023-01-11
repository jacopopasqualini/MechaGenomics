from ftplib import FTP
import subprocess

FTP_SITE = "ftp.ebi.ac.uk"
FTP_GENOMES = "/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/species_catalogue"

LOCAL_PRT_DIR = './species_proteomes/'

ftp = FTP(FTP_SITE)
ftp.login()
ftp.cwd(FTP_GENOMES)

MGYG_group = []
ftp.retrlines('NLST', MGYG_group.append)

for m in MGYG_group:
    
    FTP_MGYG = FTP_GENOMES+'/'+m
    ftp.cwd(FTP_MGYG)
    
    MGYG_group_genomes = []
    ftp.retrlines('NLST', MGYG_group_genomes.append)
    
    for g in MGYG_group_genomes: 
        
        FTP_FILE = FTP_SITE+FTP_MGYG+'/'+g+'/genome/'+g+'.faa'
        subprocess.run(['wget',FTP_FILE,'--directory-prefix',LOCAL_PRT_DIR])
                                                                    
