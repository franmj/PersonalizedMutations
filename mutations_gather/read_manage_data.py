'''
Created on May 4, 2016

@author: fran
'''
import pandas as pd
import re
import os
import subprocess
# Given a vcf file outputs a dataframe with each row representing a sample ID and every column the number of mutations in the sample within that 
# trinucleotide context

def run_command(command):
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,shell=True)
    return iter(p.stdout.readline, b'')
# Reads the vcf file and outputs a dataframe with every row representing a field in the vcf sample ID and columns with sample.id, chr, pos, ref, alt 
def read_vcf(vcf_file, sample_id):
    ##CHROM  POS     ID    REF     ALT
    f = open(vcf_file)
    table = []
    for line in  f: 
        line = line.rstrip()
        data = line.split('\t')
        if(re.match("^#", line ) != None):
            continue 
        if(len(data) ==5 ):
            if(data[3] in ("A","C","T","G") and data[4] in ("A","C","T","G")):
                table.append([sample_id,data[0],data[1],data[3],data[4]])
        if(len(data) == 4 ): 
            if(data[2] in ("A","C","T","G") and data[3] in ("A","C","T","G")):
                table.append([sample_id,data[0],data[1],data[2],data[3]])
    f.close()
    headers = names=["SAMPLE","CHROM","POS","REF","ALT"]
    df = pd.DataFrame(table, columns=headers)
   
    #table = pd.read_csv(vcf_file,sep='\t',usecols=["CHROM","POS","REF","ALT"],names=[CHROM","POS","REF","ALT)
    return df

# Given an input df with the mutations, outputs a dataframe with the number of mutations of each of the triplets
def get_triplet(df_input,sample_id,file_triplets,path_genome):
    df_triplets = create_df_triplets(file_triplets,sample_id)
    # TO DO
    df_input = df_input.sort("CHROM")
    for index, row in df_input.iterrows():
        chrom = row["CHROM"]
        pos = row["POS"]
        ref = row["REF"]
        alt = row["ALT"]
        
        fa_name =  path_genome  + chrom +".fa"
        if(not(os.path.exists(fa_name))):
            continue
        bef_sel_pos = int(pos) - 1 #base 5'
        af_sel_pos = int(pos) + 1 #base 3'
        command = "grep -v '>' %s | tr -d '\n'| cut -b%i-%i" %(fa_name,bef_sel_pos,af_sel_pos)
        #os.system(command)
        #command = 'mysqladmin create test -uroot -pmysqladmin12'.split()
        for line in run_command(command):
            print ref
            print(line)
 
    
    return df_triplets    




# create an empty dataframe with the the sample_id and the different trinucleotide mutations
def create_df_triplets(file_triplets,sample_id):
    f = open(file_triplets)
    l = [sample_id]
    values = [sample_id]
    for line in f: 
        line = line.rstrip()
        l.append(line)
        values.append(0)
        
    f.close()
    return pd.DataFrame([values],columns=l)




trinu = "/home/fran/Documents/PersonalizedProbabilities/trinucleotides.mut"    
file="/home/fran/Documents/PersonalizedProbabilities/data/294.som2.vcf"
sample="prueba"
path_genome="/home/fran/Documents/PersonalizedProbabilities/reference_genome/chromFa/"
get_triplet(read_vcf(file,sample),sample,trinu,path_genome)


