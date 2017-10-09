#!/usr/bin/env python
# coding=utf-8

import sys
import csv
#requires python 3
#<parsed.verbose> <pfam domain of interest> <e-value 0.001> 

def findDomain(pfam,dquery,evalue_input):
    for domain in pfam.split('~'):
        dtype,ddata=domain.split('(')
        evalue2=(ddata.split(" ")[2])
        evalue=float(evalue2.split('=')[1].strip(')')) #no need for float as python 3 ? [:-1] remove last charachter 
        if dquery in dtype:
            if evalue < float(evalue_input):
                return "TRUE"
#        else:
#            return "False"
    #return "False"

with open(sys.argv[1]) as fi:
    for row in csv.reader(fi, delimiter='\t'):
        if findDomain(row[1],sys.argv[2],sys.argv[3]) =="TRUE":
            print (row)    
#        if findDomain(row[1],sys.argv[2],sys.argv[3]) =="False": #Debugging
#            print ("False")
#        else: #Debugging
#            print ("Error")
