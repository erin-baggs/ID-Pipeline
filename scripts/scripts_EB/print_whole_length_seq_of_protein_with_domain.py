#!/usr/bin/env python
import sys
import csv

from ktoolu_io import readFasta

# GSMUA_Achr11P25110_001  NB-ARC(start=162, stop=447, evalue=3.4e-75)~LRR_8(start=549, stop=600, evalue=0.00014)~LRR_4(start=770, stop=802, evalue=0.12)
#If two domains in same sequence only one Fa sequence 
domain_info = dict()
valid_domains = {'DUF640'} #Python 3 {, ,} comma sep = set and {:} colon separated = dict 
with open(sys.argv[1]) as fi:
    for row in csv.reader(fi, delimiter='\t'): #changed from stdin to sys.argv
        if row:
            domains = list()
            #row[0]	GSMUA_Achr11P25110_001  
            #row[1]		NB-ARC(start=162, stop=447, evalue=3.4e-75)~LRR_8(start=549, stop=600, evalue=0.00014)~LRR_4(start=770, stop=802, evalue=0.12)
            for domain in row[1].split('~'):
                #domain		NB-ARC(start=162, stop=447, evalue=3.4e-75)
                name, coord = domain.split('(')
                if name in valid_domains:
                    domains.append({'name': name})
                    #print (name)
                    #coord 		start=162, stop=447, evalue=3.4e-75)
                    for item in coord.split(', ')[:-1]:
                        #['start=162' 'stop=447']
                        key, value = item.split('=')
                        domains[-1][key] = int(value)
            if domains:
                domain_info[row[0]] = domains


for _id, _seq in readFasta(sys.argv[2]):            
   # print(_id)    
    _id=_id.split()[0]
    if _id[1:] in domain_info : #[1:] to remove > from fasta header in readFasta tool. sequences asks if its in dictionary general
            #print(_id)
            for domain in domain_info[_id[1:]]: 
                did = '>{}:{}:{}-{}'.format(_id[1:], domain['name'], domain['start'], domain['stop'])
                #if domain['start'] - 20 > 0 and domain['stop'] + 20 < len(_seq): 
                    #dseq = _seq[domain['start'] - 20: domain['stop'] + 20]
                #if domain['start'] < 0 and domain['stop'] + 20 < len(_seq):
                    #dseq = _seq[ 0 : domain['stop'] + 20]
                #if domain['start'] > 0 and domain['stop'] + 20 > len(_seq):
                    #dseq = _seq[ domain['start'] : len(_seq)]
                #else : 
                    #dseq = _seq[ 0 : len(_seq) ]
                print(did + '\n' + _seq )
                
            #info=row[1-].split(~)
            #domain,length=info.split('(')                    

         
