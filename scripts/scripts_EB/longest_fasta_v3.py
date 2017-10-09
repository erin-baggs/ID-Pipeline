#!/usr/bin/env python
import sys

"""
gene['AT4G']={}
>>> gene['AT4G']['AT4G.2']=13
"""
#Updated: Compatible for sequence IDs with multiple dots as long as last dot seperates transcript number
#Updated: Compatible with proteomes with no dots and dots in seq ids 

gene = {}
seen = set() 
from ktoolu_io import readFasta
for _id, _seq in readFasta(sys.argv[1]):
    pid = _id.split()[0]
    #print(pid) #debug
    if  (pid.rfind('.')) > 0: 
        gid = '.'.join(pid.split('.')[:-1])
    else : 
        gid = pid  
    
    #print(gid) #debug
    if gid not in gene:
        gene[gid]={}
        gene[gid] = [pid, len(_seq), _seq]
    if float(gene[gid][1]) < len(_seq):
        #print (float(gene[gid][1]))
        #print (_seq)
        gene[gid] = [pid, len(_seq), _seq]
        #print(gid)
        #print (gene[gid]) #debug

for _id, _seq in readFasta(sys.argv[1]):
    pid = _id.split()[0]
    #print(pid) #debug
    if  (pid.rfind('.')) > 0:
        gid = '.'.join(pid.split('.')[:-1])
    if  (pid.rfind('.')) == -1 :
        gid = pid
    
    if gid not in seen:
        seen.add(gid) 
        if gid in gene: 
            print(gene[gid][0])
            print(gene[gid][2])
