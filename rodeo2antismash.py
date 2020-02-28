import os
import csv
import re
import multiprocessing as mp
from Bio import SeqIO

def rodeo_output_iterator(rod_dir, rod_dir_type):
    
    if not rod_dir_type:
        if 'main_co_occur.csv' in os.listdir(rod_dir):
            rod_dir_type = 'RODEO'
        else:
            rod_dir_type = 'RIPPER'
    
    if rod_dir_type == 'RODEO':
        # It woulod be easier to use Pandas, but this way is more memory-efficient.
        # It may be crucial when working with large datasets
        with open('%s/main_co_occur.csv' % rod_dir) as infile:
            infile.next()
            prev_seed = None
            table = []
            for row in csv.reader(infile):
                
                if prev_seed is None:
                    prev_seed = row[0]
                
                if row[0] == prev_seed:
                    table.append(row)
                
                else:
                    yield table
                    prev_seed = row[0]
                    table = [row]
                
            yield table
            
    elif rod_dir_type == 'RIPPER':
        for folder in os.listdir(rod_dir):
            if 'main_co_occur.csv' in os.listdir('%s/%s' % (rod_dir, folder)):
                with open('%s/%s/main_co_occur.csv' % (rod_dir, folder)) as infile:
                    
                    infile.next()
                    yield list(csv.reader(infile))



def rodeo_output_proccessing(table, bg_domains, n):
    
    operon_border_accs = []
    biosynthetic_genes = []
    
    
    operon_buffer = []
    prev_end = 0
    prev_strand = ''
    seed = table[0][0] # acsession of rodeo query

    for row in table:
        start = min(int(row[4]), int(row[5]))
        if (row[6] == prev_strand) and (start - prev_end < n):
            operon_buffer.append(row[3])
                    
        else:
            if seed in operon_buffer:
                operon_border_accs = [operon_buffer[0], operon_buffer[-1]]
                
            operon_buffer = [row[3]]
                    
        prev_end = max(int(row[4]), int(row[5]))
        prev_strand = row[6]
        
        if len(row) > 7:
            if row[7] in bg_domains:
                biosynthetic_genes.append(row[3])
        
        
    if seed in operon_buffer:
        operon_border_accs = [operon_buffer[0], operon_buffer[-1]]
        
    return (operon_border_accs, biosynthetic_genes)