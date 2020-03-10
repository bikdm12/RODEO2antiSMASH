import os
import csv
import re
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

class RodeoOutput(object):
    
    def __init__(self, table):
        
        self.query = table[0][0]
        self.table = table
    
    def table_proccessing(self, bg_domains, n):
    
        self.biosynthetic_genes = []
                
        operon_buffer = []
        operon_domain_buffer = []
        
        prev_end = 0
        prev_strand = ''
        
        for row in self.table:

            if len(row) > 7:
                domain = row[7]
                if row[7] in bg_domains:
                    self.biosynthetic_genes.append(row[3])
            else:
                domain = 'no_match'

            start = min(int(row[4]), int(row[5]))
            if (row[6] == prev_strand) and (start - prev_end < n):
                operon_buffer.append(row[3])
                operon_domain_buffer.append(domain)
            else:
                if self.query in operon_buffer:
                    self.operon_accs = operon_buffer
                    self.operon_domains = operon_domain_buffer
                
                operon_buffer = [row[3]]
                operon_domain_buffer = [domain]        
            
            prev_end = max(int(row[4]), int(row[5]))
            prev_strand = row[6]
                     
        if self.query in operon_buffer:
            self.operon_accs = operon_buffer
            self.operon_domains = operon_domain_buffer



def rodeo_output_iterator(rod_dir, rod_dir_type):
    
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
                    yield RodeoOutput(table)
                    prev_seed = row[0]
                    table = [row]
                
            yield RodeoOutput(table)
            
    elif rod_dir_type == 'RIPPER':
        for folder in os.listdir(rod_dir):
            if 'main_co_occur.csv' in os.listdir('%s/%s' % (rod_dir, folder)):
                with open('%s/%s/main_co_occur.csv' % (rod_dir, folder)) as infile:
                    infile.next()
                    infile_as_list = list(csv.reader(infile))
                    if len(infile_as_list) != 0:
                        yield RodeoOutput(infile_as_list)



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



def convert_gbk(gb_dir, gb_out_dir, rodeo_output, bg_domains, n, product_class):
    
    rodeo_output.table_proccessing(bg_domains, n)
    operon_border_accs = (rodeo_output.operon_accs[0], rodeo_output.operon_accs[-1])
    biosynthetic_genes = rodeo_output.biosynthetic_genes
    
    contig_edge = False
    prot_id = rodeo_output.query
    try:    
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
            cluster_qualifiers = OrderedDict([('contig_edge', str(contig_edge)), ('product', product_class)])
            cluster = SeqFeature(location = cluster_location, type = 'cluster', qualifiers = cluster_qualifiers)
            record.features = [cluster] + record.features
            
            SeqIO.write(record, '%s%s.gbk' % (gb_out_dir, prot_id), 'genbank')
    
    except Exception as e:
        print e