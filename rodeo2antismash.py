import os
import csv
import re
import multiprocessing as mp
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation



def rodeo_output_iterator(rod_dir, rod_dir_type):
    
    if rod_dir_type is None:
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
    rodeo_query = table[0][0] # acsession of rodeo query

    for row in table:
        start = min(int(row[4]), int(row[5]))
        if (row[6] == prev_strand) and (start - prev_end < n):
            operon_buffer.append(row[3])
                    
        else:
            if rodeo_query in operon_buffer:
                operon_border_accs = [operon_buffer[0], operon_buffer[-1]]
                
            operon_buffer = [row[3]]
                    
        prev_end = max(int(row[4]), int(row[5]))
        prev_strand = row[6]
        
        if len(row) > 7:
            if row[7] in bg_domains:
                biosynthetic_genes.append(row[3])
        
        
    if rodeo_query in operon_buffer:
        operon_border_accs = [operon_buffer[0], operon_buffer[-1]]
        
    return (operon_border_accs, biosynthetic_genes)



def check_if_border(feature, operon_borders):
    
    prot_id_regexp = re.compile('[A-Z]{2}_[0-9]+\.[0-9]')
        
    start_id, end_id = operon_borders
    
    if 'protein_id' in feature.qualifiers:
                    
        if feature.qualifiers['protein_id'][0]  == start_id:
            return ('start', feature.location.start + 1) # genbank is 1-based, python is 0-based
                           
        if feature.qualifiers['protein_id'][0]  == end_id:
            return ('end', int(feature.location.end))
                
    elif 'pseudo' in feature.qualifiers:        
        if 'inference' in feature.qualifiers:
                        
            inference_prot_id_search = prot_id_regexp.search(feature.qualifiers['inference'][0])
                        
            if inference_prot_id_search is not None:
                            
                inference_prot_id = inference_prot_id_search.group(0)
                            
                if inference_prot_id  == start_id:
                    return('start', feature.location.start + 1)
                        
                if inference_prot_id  == end_id:
                    return ('end', int(feature.location.end))
                    
                else:
                    if feature.qualifiers['locus_tag'][0]  == start_id:
                        return ('start', feature.location.start + 1)

                    if feature.qualifiers['locus_tag'][0] == end_id:
                        return ('end', int(feature.location.end))
                    
    return None



def convert_gbk(gb_dir, gb_out_dir, table, bg_domains, n, product_class):
    
    operon_border_accs, biosynthetic_genes = rodeo_output_proccessing(table, bg_domains, n)
    
    contig_edge = False
    prot_id = table[0][0]
        
    genbank = SeqIO.parse('%s%s.gbk' % (gb_dir, prot_id), 'genbank')
    for record in genbank: # Every file is expected to contain only one record
        
        cluster_coords = OrderedDict([('start', 1), ('end', len(record))])
        
        for feature in record.features:
            if feature.type == 'CDS':
                
                border_check = check_if_border(feature, operon_border_accs)
                if border_check is not None:
                    cluster_coords[border_check[0]] = border_check[1]

                if 'protein_id' in feature.qualifiers:
                    if feature.qualifiers['protein_id'][0] in biosynthetic_genes:
                        feature.qualifiers['sec_met'] = ['Kind: biosynthetic']
        
        start, end = cluster_coords.values()
        cluster_location = FeatureLocation(start, end)
        cluster_qualifiers = OrderedDict([('contig_edge', str(contig_edge)), ('product', 'product_class')])
        cluster = SeqFeature(location = cluster_location, type = 'cluster', qualifiers = cluster_qualifiers)
        record.features = [cluster] + record.features
        
        SeqIO.write(record, '%s%s.gbk' % (gb_out_dir, prot_id), 'genbank')