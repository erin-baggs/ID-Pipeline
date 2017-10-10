#!/usr/bin/env python
from __future__ import print_function
import sys
import csv
import re
import os
import subprocess 
from blast import *
from ktoolu_io import readFasta


# GSMUA_Achr11P25110_001  NB-ARC(start=162, stop=447, evalue=3.4e-75)~LRR_8(start=549, stop=600, evalue=0.00014)~LRR_4(start=770, stop=802, evalue=0.12)
#If two domains in same sequence only one Fa sequence 
gid_nlrid = list()
id_names = set()
domain_info = dict()
nlrids = set()
valid_domains = {'NB-ARC','LRR','NB-LRR','TIR','RPW8','DUF3542', 'TIR_2', 'DUF640'} #Add other DUFs I have charachterised  #Python 3 {, ,} comma sep = set and {:} colon separated = dict 
with open(sys.argv[1]) as fi:
    for row in csv.reader(fi, delimiter='\t'): #changed from stdin to sys.argv
        if row:
            gid_nlrid.append({'nlrid': row[0]}) #List of gids associated with NLR-IDs   
            domains = list()
            #row[0]	GSMUA_Achr11P25110_001  
            #row[1]		NB-ARC(start=162, stop=447, evalue=3.4e-75)~LRR_8(start=549, stop=600, evalue=0.00014)~LRR_4(start=770, stop=802, evalue=0.12)
            for domain in row[1].split('~'):
                #domain		NB-ARC(start=162, stop=447, evalue=3.4e-75)
                name, coord = domain.split('(')
                if name not in valid_domains and not (name.startswith('LRR') or  name.startswith('LRV')):
                    domains.append({'name': name})
                    id_names.add(name)
                    #print (name)
                    #coord 		start=162, stop=447, evalue=3.4e-75)
                    for item in coord.split(', ')[:-1]:
                        #['start=162' 'stop=447']
                        key, value = item.split('=')
                        domains[-1][key] = int(value)
            if domains:
                domain_info[row[0]] = domains
                nlrids.add(row[0])

#Need to keep list of IDs' with NLR-IDs
#Now have list of domains - need to extract all proteins with these domains from parsed.verbose

with open(sys.argv[2]) as fi2:
    for row in csv.reader(fi2, delimiter='\t'): #changed from stdin to sys.argv
        if row:
            domains = list()
            for domain in row[1].split('~'):
                name, coord = domain.split('(')
                if name in id_names:
                    domains.append({'name': name})
                    #coord              start=162, stop=447, evalue=3.4e-75)
                    for item in coord.split(', '):
                        #['start=162' 'stop=447']
                        key, value = item.strip(')').split('=')
                        cast = float if key == 'evalue' else int 
                        domains[-1][key] = cast(value)
                    if domains[-1]['evalue'] > 0.001 : #Filter proteome annotations to be >0.001
                        del domains[-1]
            if domains:
                domain_info[row[0]] = domains


seen_id={}
for _id, _seq in readFasta(sys.argv[3], headless=True):            
    # print(_id)    
    _id=_id.split()[0]
    #counter=0
    gid ='initial'
    pid = 'initial'
    # bait_data['intial']=['intial',0]
    #print (_id)
    if _id in domain_info : #[1:] to remove > from fasta header in readFasta tool. sequences asks if its in dictionary general
        #Here I could filter for 1 per primary transcript
        for domain in domain_info[_id]:
            gid = ".".join(_id.split('.')[:-1])
            if gid not in seen_id:
                seen_id[gid] = [_id, len(_seq)]
            elif len(_seq) > seen_id[gid][1]:
                seen_id[gid] = [_id, len(_seq)]
#print (seen_id)
#print (domain_info)
with open('query.fa', 'w') as query_out, open('db.fa', 'w') as db_out :
    files = dict()
    for _id, _seq in readFasta(sys.argv[3], headless=True):
       # print(_id) 
        _id = _id.split()[0]
        gid = ".".join(_id.split('.')[:-1])
        if _id in domain_info and gid in seen_id and seen_id[gid][0] == _id : #[1:] to remove > from fasta header in readFasta tool. sequences asks if its in dictionary general
            print(_id)
            #Here I could print the indvidual domain/NLR-ID fasta files  
            for domain in domain_info[_id]:
                if files.get(domain['name'], None) is None:
                    files[domain['name']] = domain['name'] +'_domain.fa', open(domain['name'] +'_domain.fa', 'w')
                isID, out = "_DOMAIN", db_out
                if _id in nlrids:
                    isID, out = "_NLR-ID", query_out
                did = '>{}~{}_{}-{}{}'.format(_id, domain['name'], domain['start'], domain['stop'], isID)
                if domain['start'] - 20 > 0 and domain['stop'] + 20 < len(_seq): 
                   dseq = _seq[domain['start'] - 10: domain['stop'] + 10]
                if domain['start'] < 0 and domain['stop'] + 10 < len(_seq):
                   dseq = _seq[ 0 : domain['stop'] + 10] 
                if domain['start'] > 0 and domain['stop'] + 10 > len(_seq):
                   dseq = _seq[ domain['start'] : len(_seq)]
                else : 
                   dseq = _seq[ 0 : len(_seq) ]
                print(did + '\n' + dseq , file = out )
                print(did + '\n' + dseq , file = files[domain['name']][1])
             
    for k in files:
        try:
            files[k][1].close()
        except: 
            pass


blast_seen=set()
filtered_qid_sid=dict()
with open('raw.blast', 'w') as blast_out, open('filtered.blast', 'w') as blast_filtered:

    stdout, stderr = makeBlastdb('db.fa')
    if not stderr :
        print ('working')
        for HSP in runBlast('query.fa', BLAST_CMD, 'db.fa', 'blastp', nthreads=1): #default nthreads = 8  
            print(*HSP, sep = '\t', file = blast_out)
            if HSP.query not in blast_seen:
            #if float(HSP.evalue) < 1e-19:
                blast_seen.add(HSP.query)
                print(*HSP, sep = '\t', file = blast_filtered)
                subject = '.'.join(HSP.subject.split('~')[0].split('.')[:-1])
                query = '.'.join(HSP.query.split('~')[0].split('.')[:-1])
                filtered_qid_sid[subject] = filtered_qid_sid.get(subject, list())#new
                filtered_qid_sid[subject].append(query) 
    else: #old
        print(stderr)

print (filtered_qid_sid)
        
GFF_dict= {}
with open(sys.argv[4]) as GFF, open('chromo_cord.gff', 'w') as chromo_coord:
    for row in csv.reader(GFF, delimiter='\t'): 
        if row[2] == 'gene':
            #print (row)
            #Chr1    TAIR10  gene    3631    5899    .       +       .       ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
            for item in row[8].split(';'):
                key, value = item.split('=')
                GFF_dict[key] = str(value)
            if GFF_dict['ID'] in filtered_qid_sid : #regex as GFF maybe different #Can I cat together GFFs
                for query in filtered_qid_sid[GFF_dict['ID']]:
                    print (*(row + [query]), sep = '\t' , file = chromo_coord)

 
"""
source t_coffee-9.03.r1318

QUERY=$1
OUT=$2

t_coffee ../../data/*_domain.fa -mode mcoffee -outfile ../../data/AIG_1.mcoffee.fa -output fasta_aln
"""
done_files=list()
TCOFFEE_CMD = 't_coffee {} -mode mcoffee -outfile {} -output fasta_aln; touch {};'
TCOFFEE_SBATCH = 'sbatch -p {} -c {} --mem {} --wrap "{}" -J EB_Pipe_Tcoffee'
for domain_fasta in files: 
    tcoffee_cmd = TCOFFEE_CMD.format(files[domain_fasta][0], files[domain_fasta][0]+'.aln', files[domain_fasta][0]+'tcoffee.done')  
    sbatch_cmd = TCOFFEE_SBATCH.format('ei-medium', 1, '2GB', tcoffee_cmd)
    done_files.append(files[domain_fasta][0]+'tcoffee.done')
    pr = subprocess.Popen(sbatch_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    out, err = pr.communicate()
    print(out, err, sep = '\n')
