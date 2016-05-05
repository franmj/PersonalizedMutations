from __future__ import with_statement
import os
from string import join
import glob
import itertools
from count import *
import pandas as pd
from pandas import Series

path = '/scratch/Samples/inputData/'
names=[]
dict_input = {}
list_names=[]
reference_list = []

#extract the name of the sample, the mutated base, and neighbors of that base, output a matrix
for sname in os.listdir(path):
    if sname.endswith("43.som2.vcf"):
        shandler = open(path + sname)
        key = sname
        while True:
            try:
                line = shandler.next()
            except StopIteration:
                break
            if line[0].islower():
                new_line = line.strip().split('\t') #separate elements line '\t'
                chromosome = new_line[0] #select chromosome
                position = new_line[1] #position of mutated base
                ref_alt = new_line[2:] #base of reference and altered
                ref = ref_alt[0]
                alt = ref_alt[1]
                if len(ref) == 1 and len(alt) == 1 or "," in alt: #only take into account substitutions, avoid #indels
                    mut_type = join(ref + alt,">")
                    if ',' in mut_type:
                        alt1 = alt[0]
                        alt2  = alt[2]
                        mut_type = '>'.join(ref + alt1) + ":" + alt2
                    if key in dict_input: #dict_input: key = sample name; values = chr, position, ref, alt.
                        dict_input[key].append((chromosome,position,mut_type))
                    else:
                        dict_input[key] = []
                        dict_input[key].append((chromosome,position,mut_type))
        shandler.close()

#look for the neighbor bases of the altered
path_genome = '/scratch/chromFa/*.fa'
dict_final ={}
#open chromosome file depending on name, when opened look for position mutated base
for key,value in dict_input.items():
    sel_chromosome = [chrom[0] for chrom in value] #select chromosome where mutated base is found
    sel_pos = [pos[1] for pos in value] #select position to once opened file look neighbor bases
    sel_mut_type = [pos[2] for pos in value] #select mutation type
    for schromosome,spos,smuttype in itertools.izip(sel_chromosome,sel_pos,sel_mut_type): #loop chromosome and position
        for filename in glob.glob(path_genome): #list all the file names in directory
            files = os.path.basename(filename)[:-3] #from path keep only name of the chr to match
            if schromosome == files:
                ghandler = open(filename) #open the file and the name chr
                while True:
                    try:
                        bef_sel_pos = int(spos) - 1 #base 5'
                        af_sel_pos = int(spos) + 1 #base 3'
                        command = "grep -v '>' %s | tr -d '\n'| cut -b%i-%i > subtype3.txt" %(filename,bef_sel_pos,af_sel_pos)
                        p = os.system(command)
                        file = open("subtype3.txt")
                        for output in file:
                            output = output[:-1]
                            if key in dict_final:
                                dict_final[key].append((smuttype,output))
                            else:
                                dict_final[key] = []
                                dict_final[key].append((smuttype,output))
                        break
                    except StopIteration:
                        break
                ghandler.close()

#
# # array,sample_names = count(dict_final)
#
# # scipy.io.savemat('/home/sgalan/Desktop/all_CLL_WGS_CNAG.mat', mdict={'originalGenomes':array})
# # scipy.io.savemat('/home/sgalan/Desktop/all_CLL_WGS_CNAG_sample_names.mat', mdict={'Names':sample_names})
#
