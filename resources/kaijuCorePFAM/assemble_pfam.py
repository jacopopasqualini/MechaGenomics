import pandas as pd
import os

#import matplotlib
#import matplotlib.pyplot as plt

LOCAL_DIR = './'
PROTEOMES_DIR = os.path.join(LOCAL_DIR,'proteomes')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# build pfam abundance matrix

missed = []

Abundances = pd.DataFrame()

for p in os.listdir(PROTEOMES_DIR):
        
    print(' > > > ',p)
    
    proteome_f = os.path.join(PROTEOMES_DIR,p)

    if os.path.exists( proteome_f ):
    
        try:
            proteome = pd.read_csv(proteome_f,compression='gzip',skiprows=3, header=None, sep='\t')
            proteome = proteome[5].value_counts()
            proteome = pd.Series(proteome,name=p.replace('.tsv',''))
            Abundances = pd.concat([Abundances,proteome],axis=1)
        
        except EOFError:
            print('EOFError: ',p)
            missed.append(p)
            
    else:
        print('Missed file: ',p)
        missed.append(p)

Abundances.index=Abundances.index.rename('PF')

Abundances.to_csv('PFAM_proteome_matrix.csv',index=True)

file=open('./sheets/missed_proteomes.txt','w')
for items in missed: file.writelines(items+'\n')
file.close()
