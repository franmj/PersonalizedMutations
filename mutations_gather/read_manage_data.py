'''
Created on May 4, 2016

@author: fran
'''
import pandas as pd
import re
import os
import subprocess
import glob
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
    df.set_index('SAMPLE')
   
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
    l = ['SAMPLE']
    values = [sample_id]
    for line in f: 
        line = line.rstrip()
        l.append(line)
        values.append(0)
        
    f.close()
    df = pd.DataFrame([values],columns=l)
    df.set_index('SAMPLE')
    return df

#Directly loads the trinucleotide frequency of a file inot an empty dataframe
def load_counts(file,sample_id):
    
    f=open(file)
    col = ['SAMPLE']
    values = [sample_id]
    for line in f: 
        line = line.rstrip()
        data = line.split('\t')
        if(len(data[0]) == 7):
            col.append(str(data[0][0]+data[0][2]+data[0][6]))
        else:
            col.append(data[0])
        values.append(int(data[1]))
 
    f.close()
    df = pd.DataFrame([values],columns=col)
    df.set_index('SAMPLE')
    return df

#create a file with the number of trinucleotide in the genome/exome. Very time consuming. 

def create_file_frecuencies(list_trinucleotides,path_fastas,file_output):
    
    map_values = {}
    for colname in list_trinucleotides:
        if(len(colname)!=7):
            continue
        colname = colname[0]+colname[2]+colname[6]
        value = 0
        print path_fastas+'/chr*.fa'
        for filename in glob.iglob(path_fastas+'/chr*.fa'):
            command = "grep -io "+colname+" "+filename+" | wc -l > " + path_fastas + "/tmp"
            print command
            os.system(command)
            f= open(path_fastas + "/tmp")
            for line in f:
                line = line.rstrip()
                v = int(line)
            f.close()
            value = int(v) + value
            print str(value)
            
        map_values[colname] = value
    #for file,dir,subdir
    f = open(file_output,'w')
    for key in map_values.keys():
        f.write(key + '\t' + str(map_values[key]) + "\n")
        
        
    
    
    
    f.close()



# Get a dataframe with the number of counts per trinucleotide in the exome or genome
def get_count_tri(type,sample_id,file_triplets):
    
    if("exome" in type):
        #df = create_df_triplets_basic(file_triplets,sample_id)
        df_norm = load_counts(file_exomes,sample_id)
    elif("genome" in type):
        #df = create_df_triplets_basic(file_triplets,sample_id)
        #print df
        df_norm = load_counts(file_genomes,sample_id)
    else: 
        return False
    return df_norm


file_exomes = "/Users/fran/Documents/Work/PersonalizedMutations/data/test.out"
file_genomes = "/Users/fran/Documents/Work/PersonalizedMutations/data/test.out"
trinu = "/home/fran/Documents/PersonalizedProbabilities/trinucleotides.mut"    
trinu = "/Users/fran/Documents/Work/PersonalizedMutations/trinucleotides.mut"
file="/home/fran/Documents/PersonalizedProbabilities/data/294.som2.vcf"
sample="prueba"
type = "genome"
path_genome="/home/fran/Documents/PersonalizedProbabilities/reference_genome/chromFa/"
def get_dataframe_mutations(sample,file,type="genome"):
  
    df_mutations = get_triplet(read_vcf(file,sample),sample,trinu,path_genome)
    df_genome = get_count_tri(type,sample,trinu)
    return df_mutations,df_genome


df_mut,df_genome = get_dataframe_mutations(sample,"/Users/fran//Documents/Work/PersonalizedMutations/data/294.som2.vcf")
print df_mut
print df_genome
#columns = create_df_triplets("/Users/fran/Documents/Work/PersonalizedMutations/trinucleotides.mut",sample).columns

#create_file_frecuencies(columns,"/Users/fran/Documents/Work/PersonalizedMutations/","/Users/fran/Documents/Work/PersonalizedMutations/data/test.out")