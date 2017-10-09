
#!/usr/bin/env python
# coding=utf-8

import sys
import csv
import re


"""
# BLASTN 2.2.28+
# Query: TIL0100004529-RA_subseq_538-940 coordinates 538 - 940 of TIL0100004529-RA
# Database: ../current_bait_set/baits_high_HSP_removed.fa
# Fields: query id, subject id, % identity, alignment length, mismatches, q. start, q. end, s. start, s. end, evalue, bit score, query length, subject length, identical, positives, gap opens, gaps, % positives, % subject coverage, % hsp coverage, subject strand, query seq, subject seq
# 1 hits found
TIL0100004529-RA_subseq_538-940 lNl|Chr4_nlr_1_1122135-1126639:+_69417  100.00  120     0       49      168     1       120     2e-58    222    403     120     120     120     0       0
       100.00  30      30      plus    GTTACCGATGGTGATCTTATGGTCGTTCCCATTGTTGGAATGGGAGGACTAGGCAAGACCACACTCGCTCAGCTCATCTACAACGACCCTCAGGTCAAGGAGCACTTCCATCTACTAAAG        GTTACCGATGGTGATCTTATGGTCGTTCCCATTGTTGGAATGGGAGGACTAGGCAAGACCACACTCGCTCAGCTCATCTACAACGACCCTCAGGTCAAGGAGCACTTCCATCTACTAAAG

TIL0100046533-RA_subseq_151-325 NLR
TIL0100049719-RA_subseq_466-676 NLR
TIL0100052409-RA_subseq_544-727 NLR
"""
counter=0
with open(sys.argv[1]) as fi:
    seen_id={}
    #counter=0
    gid ='initial'
    pid = 'initial'
    # bait_data['intial']=['intial',0]
    for row in csv.reader(fi, delimiter='\t'):
        pid = str(row[0])
        gid = str(row[0].rsplit('.', 1)[0])
        if gid not in seen_id:
            seen_id[gid] = [pid]
                       
for k in seen_id:
    output=seen_id[k]
    output=str(output)
    #for i in output:
    # print('\t'.join(i))
    #print (output)
    print (seen_id[k])

                
#with open(sys.argv[2]) as fi1:
    #for row in csv.reader(fi1, delimiter=','):
        #if row[0] in bait_data:
           # printing_string=[]
           # printing_string=bait_data[row[0]]
            #printing_string=str(printing_string)
            #try:
                #print(printing_string+'\t'+row[1]+'\t'+row[2])
            #except:
                #print(printing_string+'\t'+row[1])
            #print ("YES")
        #else:
            #print ("ERROR")
            


"""
with open(sys.argv[3]) as fi2:
    qid='intial'
    for row in csv.reader(fi2, delimiter='\t'):
        if row[0].startswith('# BLASTN 2.6.0+'):
            bait_data[qid].append([counter, 'Sitalica'])
            counter = 0
        else:
            if row[0].startswith('# Query'):
                qid=row[0].split('|')[1]
                #if row[0].startswith('# 0 hits found'):
                # bait_data[qid].append([ '0', 'Sitalica'])
                # if re.match('#.{1,3} hits found', row[0]):
                #hits=row[0].split(' ')[1]
            elif not row[0].startswith('#'):
                q_len, aln_len, sub_len, pid, gaps, mismatch = map(float, ( row[11], row[3], row[12], row[2], row[16] , row[4]))
                align_value=(120-(gaps + mismatch))/120
                if align_value >=0.79:
                   counter = counter+1
                   #print (counter)
                   #bait_data[qid].append([hits, 'Sitalica'])
                   #if qid not in bait_data:
                   # bait_data[qid].append(['0', 'Sitalica'])

bait_data[qid].append([counter, 'Sitalica'])
for k in bait_data:
    output=bait_data[k]
    output=str(output)
    #for i in output:
       # print('\t'.join(i))
    #print (output)
    print (bait_data[k])
"""            
"""
                #debug.append('yes')
    #return debug
    return bait_data
     

with  open(sys.argv[1]) as fi:
    output=bait_dict(fi)

with  open(sys.argv[2]) as fi1:
    output2=bait_dict(fi1)

with  open(sys.argv[3]) as fi2:
    output3=bait_dict(fi2)


    for k in output2:
        list1=output2[k]
        #if k in output[k]:
        output[k].append[list1]
    for k in output3:
        list1=output3[k]
        output[k].append[list1]

print(output)

"""
