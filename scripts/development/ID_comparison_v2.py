#!/usr/bin/env python
from __future__ import print_function
import sys
import csv
import re
import os
import subprocess 
import time
import glob
from blast import *
from ktoolu_io import readFasta



##FUNCTION : Identify all domains matching ID
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



##FUNCTION: Filter ID matching  domains for e-value threshold 
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

##FUNCTION : Record only longest fasta transcript of ID matching domains
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


## FUNCTION: Write out ID matching domains fasta 
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

## FUNCTION : BLAST NLR-ID vs ID domain fasta 
blast_seen=set()
alignment_domains=set()
filtered_qid_sid=dict()
with open('raw.blast', 'w') as blast_out, open('filtered.blast', 'w') as blast_filtered:

    stdout, stderr = makeBlastdb('db.fa')
    if not stderr :
        print ('working')
        for HSP in runBlast('query.fa', BLAST_CMD_ALN, 'db.fa', 'blastp', nthreads=1): #default nthreads = 8  
            print(*HSP, sep = '\t', file = blast_out)
            domains_2 = HSP.subject
            domains_3 = domains_2.split('~')[1]
            domains_4 = domains_3.split('_')[0]
            alignment_domains.add(domains_4)
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
print (alignment_domains)

"""
##FUNCTION : Seqeuences to include in alignment
# HELP 1. dont think I can use query need to do seperate for each domain as above 
#  2. how would i then read input fasta for that blast and will I need to split/ join ??? 
blast_seen=set()
with open('to_aln.blast') as blast_aln: 
    for HSP in runBlast('query.fa', BLAST_CMD_ALN, 'db.fa', 'blastp', nthreads=1):
        if HSP.subject not in blast_seen_2:
            blast_seen_2.add(HSP.subject)
            #print(HSP.subject, sep = '\t', file = to_aln.blast)
        if HSP.query not in blast_seen_2:
            blast_seen_2.add(HSP.query)
            #print(HSP.query, sep = '\t', file = to_aln.blast )
 
    for _id, _seq in readFasta(, headless=True):
        if _id      
"""

## FUNCTION : Print GFF co-ordinates of paralog to ID         
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

#NEW WITHOUT CSHU 
##Function get FASTA from db.fa of all target domains
#with open(db.fa) as fasta_db:
files_2=dict()
#print(alignment_domains)
for _id, _seq in readFasta('db.fa', headless=True):
    #print(_id)
    print(alignment_domains)
   # if _id in alignment_domains:
    domain_long = _id.split('~')[1]
    domain_short = domain_long.split('_')[0]
    print(domain_short)
    if domain_short in alignment_domains:
        if files_2.get(domain_short, None) is None:
            files_2[domain_short] = domain_short + 'pre-aln.fa', open(domain_short + 'pre-aln.fa', 'w') 
        print('>'+_id + '\n' + _seq , file = files_2[domain_short][1])
        print ('db.fa seen')
## NEED TO ADD QUERY TO CORRECT OUT FILE 
for _id, _seq in readFasta('query.fa', headless=True):#Can I access newly created file this way? 
    q_domain_long = _id.split('~')[1]
    q_domain_short = q_domain_long.split('_')[0]
    if files_2.get(q_domain_short, None) is None:
        pass
    else:
        print('>'+_id + '\n' + _seq , file = files_2[q_domain_short][1])
        print ('query.fa seen')

for k in files_2:
    try:
        files_2[k][1].close()
    except:
        pass


## FUNCTION : Multiple sequence alignment of each *domain.fa to nlr_id
## NEED TO WORK ON TCOFFEE SLURM AS WRONG FORMAT !! 
start_files=list()
done_files=list()
TCOFFEE_CMD = 'touch {}; t_coffee {} -mode mcoffee -outfile {} -output fasta_aln; touch {};'
TCOFFEE_SBATCH = 'sbatch -p {} -c {} --mem {} --wrap "{}" -J EB_Pipe_Tcoffee'
list_files_2 = glob.glob("*pre-aln.fa") #glob is imported for regular expression finding. Can add path to glob if create in subdirectory sub/"." everything in here print  

for f in list_files_2:   
    tcoffee_cmd = TCOFFEE_CMD.format(f+'tcoffee.start', f, f+'.aln', f+'tcoffee.done')  
    sbatch_cmd = TCOFFEE_SBATCH.format('ei-medium', 1, '2GB', tcoffee_cmd)
    done_files.append(f+'tcoffee.done')
    start_files.append(f+'tcoffee.start')
    pr = subprocess.Popen(sbatch_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    out, err = pr.communicate()
    print(out, err, sep = '\n')

start_time = time.time()
unfinished = list(done_files)
while True:
    unfinished = [f for f in unfinished if not os.path.exists(f)]
    if not unfinished or (time.time()-start_time)/60:
        break

    time.sleep(300) 

"""    
## Function : Curation of alignment 
???


## Function : RAXML job submissions
#RAXML job code 
#Sort retrieve alignment files, file naming and time
raxml_start_file = list()
raxml_done_file = list()
RAXML_CMD = 'touch {}; raxmlHPC-MPI-SSE3 -f a -x 1123 -p 2341 -# 100 -m PROTCATJTT -s $1  -n {}; touch {};' #-n = output name #UNFINISHED
RAXML_SBATCH = 'sbatch -p {} -c {} --mem {} --wrap "{}" -J EB_Pipe_RAXML'
for ????
    raxml_cmd = #raxmlHPC-MPI-SSE3 -f a -x 1123 -p 2341 -# 100 -m PROTCATJTT -s $1  -n $(basename $1)raxml
    sbatch_cmd = RAXML_SBATCH.format('ei-long', 4, '48GB', raxml_cmd)
    raxml_done_file.append(?+'raxml.done')
    raxml_start_file.append(?+'raxml.start')
    pr = subprocess.Popen(sbatch_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = pr.communicate()
    print(out, err, sep = '\n')

start_time = time.time()
unfinished = list(done_files)
while True:
    unfinished = [f for f in unfinished if not os.path.exists(f)]
    if not unfinished or (time.time()-start_time)/60:
        break

    time.sleep(??days??)
"""
