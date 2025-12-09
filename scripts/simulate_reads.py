import argparse
import os
import random
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO

bases = {1:'a', 2:'g', 3:'c', 4:'t'}

def insertion(string, i):
    start = string[:i]
    end = string[i:]
    to_add_string = ""
    for i in range (random.randint(1,10)):
        to_add_string += bases[random.randint(1,4)]
    return (start + to_add_string + end), to_add_string

def deletion(string, i, length):
    start = string[:i]
    end = string[i+length:]
    return start + end

def SNP(string, i): 
    current = string[i].lower()
    done = False
    while done == False:
        
        new = bases[random.randint(1,4)]
        if new != current:
            done = True
    return (string[:i] + new + string[i+1:]), current, new, i

def check_index(string, i):
    return 0 <= i < len(string) and string[i] == string[i].upper()

def mutate_genome(genome):
    
    original = genome
    mutations = pd.DataFrame(columns=["index","original","new"])
    
    ##SNPS
    for i in range(300):
        found_good_index = False
        while not found_good_index:
            x = random.randint(0, len(genome)-1)
            if check_index(genome, x): 
                found_good_index = True
        
        genome, current, new, index = SNP(genome, x)
        mutations.loc[len(mutations)] = [index, current, new]
        
    ##InDels
    indicies = pd.DataFrame(columns=["index", "which"])
    for i in range(20):
        in_or_del = random.randint(1,2)
        
        #if insertion
        if in_or_del == 1:
            
            found_good_index = False
            while not found_good_index:
                x = random.randint(0, len(genome)-1)
                if check_index(genome, x) and x not in indicies["index"].values:
                    found_good_index = True #if insertion not on SNP
            
            indicies.loc[len(indicies)] = [x, None]
         
        #if deletion
        else:
            
            length = random.randint(1,10)
            found_good_index = False
            while not found_good_index:
                x = random.randint(0, len(genome)-length)
                found = True
                for i in range(length):
                    
                    if not check_index(genome, x+i):
                        found = False
                    
                if found == True and x not in indicies["index"].values:
                    found_good_index = True
                
            indicies.loc[len(indicies)] = [x, length]   
        
    #make changes based on choices
    indicies = indicies.sort_values(by=["index"]).reset_index(drop=True)
    indicies["new_index"] = indicies["index"].copy()

    for row in range(len(indicies)):
        l = indicies.at[row, "which"]
        
        if pd.notna(l):#deletion
            l = int(l)
            curr_i = int(indicies.at[row, "new_index"])
            genome = deletion(genome, curr_i, l)
            indicies.loc[row+1:, "new_index"] -= l
            
            o_index = int(indicies.at[row, "index"])  
            deleted_seq = original[o_index:o_index + l]
            mutations.loc[len(mutations)] = [o_index, deleted_seq, ""] 
            
        else:
            curr_i = int(indicies.at[row, "new_index"])
            genome, to_add = insertion(genome, curr_i)
            indicies.loc[row+1:, "new_index"] += len(to_add)
            
            o_index = int(indicies.at[row, "index"])  # original coordinate
            mutations.loc[len(mutations)] = [o_index, "", to_add]
    
    return mutations, genome

def parse_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference", required=True, help="Reference genome FASTA for simulation.")
    parser.add_argument("-n", "--name", required=True, help="Name for output files.")
    parser.add_argument("-o", "--outdir", required=True, help="Directory to write output files into.")

    return parser.parse_args()

def write_to_fasta(reads, f_name):
    
    output_file = open(f_name, "w")
    for i in range(len(reads)):
        output_file.write("@simulated_read_"+str(i+1)+"\n")
        output_file.write(reads[i] + "\n")
        output_file.write("+\n")
        qual = ""
        for j in range(len(reads[i])):
            qual += "I"
        output_file.write(qual+"\n")
    output_file.close()

def main():
    
    args = parse_args()
    reference = args.reference
    prefix = args.name
    outdir = args.outdir

    os.makedirs(args.outdir, exist_ok=True)
    mutations_file = os.path.join(args.outdir, f"{prefix}_mutations.csv")
    fastq1_file = os.path.join(args.outdir, f"{prefix}_fastq1.fq")
    fastq2_file = os.path.join(args.outdir, f"{prefix}_fastq2.fq")
    
    file = open(reference)
    sequence = ""
    for line in file:
        if line[0] != '>' and line!="\n":
            sequence += line.split('\n')[0]
     
    mutations, mutated_genome = mutate_genome(sequence)
    mutations.to_csv(mutations_file, index=False)
    
    seq = Seq(mutated_genome)
    reverse_comp = str(seq.reverse_complement())
    
    read_length = 100
    step = 6
    foward_reads = []
    reverse_reads = []
    reverse_comp = str(Seq(mutated_genome).reverse_complement())

    for start in range(0, len(mutated_genome) - read_length, step):
        end = start + read_length
        foward_reads.append(mutated_genome[start:end])
        reverse_reads.append(reverse_comp[start:end])
        
    write_to_fasta(foward_reads, fastq1_file)
    write_to_fasta(reverse_reads, fastq2_file)

if __name__ == "__main__":
    main()