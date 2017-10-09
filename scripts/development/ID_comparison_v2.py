#!/usr/bin/env python
import sys
import csv
import re

from ktoolu_io import readFasta


# GSMUA_Achr11P25110_001  NB-ARC(start=162, stop=447, evalue=3.4e-75)~LRR_8(start=549, stop=600, evalue=0.00014)~LRR_4(start=770, stop=802, evalue=0.12)
#If two domains in same sequence only one Fa sequence 
gid_nlrid = list()
id_names = list()
domain_info = dict()
valid_domains = {'NB-ARC','LRR','NB-LRR','TIR','RPW8','DUF3542'} #Add other DUFs I have charachterised  #Python 3 {, ,} comma sep = set and {:} colon separated = dict 
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
                    id_names.append(name)
                    #print (name)
                    #coord 		start=162, stop=447, evalue=3.4e-75)
                    for item in coord.split(', ')[:-1]:
                        #['start=162' 'stop=447']
                        key, value = item.split('=')
                        domains[-1][key] = int(value)
            if domains:
                domain_info[row[0]] = domains

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
                    for item in coord.split(', ')[:-1]:
                        #['start=162' 'stop=447']
                        key, value = item.split('=')
                        domains[-1][key] = int(value)
            if domains:
                domain_info[row[0]] = domains

for _id, _seq in readFasta(sys.argv[3]):            
   # print(_id)    
    _id=_id.split()[0]
    if _id[1:] in domain_info : #[1:] to remove > from fasta header in readFasta tool. sequences asks if its in dictionary general
            #print(_id)
            #Here I could filter for 1 per primary transcript 
            for domain in domain_info[_id[1:]]: 
                did = '>{}_{}_{}-{}'.format(_id[1:], domain['name'], domain['start'], domain['stop'])
                if domain['start'] - 20 > 0 and domain['stop'] + 20 < len(_seq): 
                    dseq = _seq[domain['start'] - 10: domain['stop'] + 10]
                if domain['start'] < 0 and domain['stop'] + 10 < len(_seq):
                    dseq = _seq[ 0 : domain['stop'] + 10]
                if domain['start'] > 0 and domain['stop'] + 10 > len(_seq):
                    dseq = _seq[ domain['start'] : len(_seq)]
                else : 
                    dseq = _seq[ 0 : len(_seq) ]
                print(did + '\n' + dseq )
                
            #info=row[1-].split(~)
            #domain,length=info.split('(')                    

         
