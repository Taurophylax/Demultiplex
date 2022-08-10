import gzip
import re
import argparse as ap
import matplotlib.pyplot as plt
import os
import itertools as it
import numpy as np

parser = ap.ArgumentParser(description='demux.py -i Read1.fq.gz Index1.fq.gz Index2.fq.gz Read2.fq.gz index_list.txt')
parser.add_argument('-i', '--infiles', help='In-files (fq.gz)', nargs=5, type=str)
args = parser.parse_args()
inFiles = args.infiles
in_R1file = inFiles[0]
in_i1file = inFiles[1]
in_i2file = inFiles[2]
in_R2file = inFiles[3]
in_indexes = inFiles[4]

indexes = [x.split('\t')[4].strip() for x in open(in_indexes).readlines()] #read index file and pull 5th column index sequences
indexes.pop(0) #remove header (line 1) of index file
indexperms = it.product(indexes,repeat=2) #create cartesian product of index pairs
indexpermcounts = {}
for a in indexperms:                #fill perm counts with 0
    a = a[0] + "-" + a[1]
    indexpermcounts[a] = 0
R1_matchfiles = {}
R2_matchfiles = {}
for i in indexes:
    R1_matchfiles[i] = open(f"OUT/R1_{i}.fq", "a")
    R2_matchfiles[i] = open(f"OUT/R2_{i}.fq", "a")
R1_hopped = open("OUT/R1_hopped.fq", "a")
R2_hopped = open("OUT/R2_hopped.fq", "a")
R1_unknown = open("OUT/R1_unknown.fq", "a")
R2_unknown = open("OUT/R2_unknown.fq", "a")

def revcomp(DNA: str) -> str:
    nt = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    cDNA = "".join(nt[i] for i in reversed(DNA))   #reverse DNA string and swap based on k-v pairs
    return cDNA

def repair(i1, i2) -> [str,str]:    #tries to repair N's by comparing with index match, cannot repair N-pairs (obviously)
    N = "N" # because N is N
    Ncount = 0
    N_i1 = [i.start() for i in re.finditer(N, i1)]  #iterate through i1 and add index (to list) where N occurs (if any)
    N_i2 = [i.start() for i in re.finditer(N, i2)]
    i1 = list(i1) #convert strings to a list for more iteration
    i2 = list(i2)
    Ncount = len(N_i1) + len(N_i2) #get total Ns for both strings
    if Ncount > 2:                       #automatically fail if there are more than 2 N's total
        return(''.join(i1), ''.join(i2)) #return i1 and i2 as they were originally and break 
    for x in N_i1:  #repair i1 with i2
        if i1[x] != i2[x]:
            i1[x] = i2[x]
    for x in N_i2:  #repair i2 with i1
        if i2[x] != i1[x]:
            i2[x] = i1[x]
    return(''.join(i1), ''.join(i2)) #return i1 and i2 as a string   

def getrecords(indexes: list):
    recordnum = 1
    R1 = gzip.open(in_R1file, "rt") #read R1
    I1 = gzip.open(in_i1file, "rt") #read I1
    I2 = gzip.open(in_i2file, "rt") #read I2
    R2 = gzip.open(in_R2file, "rt") #read R2
    statR1match, statR2match, statR1hopped, statR2hopped, statR1unknown, statR2unknown = 0,0,0,0,0,0
    while True:
        recordnum += 1
        recordR1, recordI1, recordI2, recordR2 = [], [], [], []

        for i in range(0,4):  #pull record and append each line to a temporary list
            recordR1.append(R1.readline().strip())
            recordI1.append(I1.readline().strip())
            recordI2.append(I2.readline().strip())
            recordR2.append(R2.readline().strip())

        if recordR1[0] == "":   #ARE WE DONE?
            return statR1match, statR2match, statR1hopped, statR2hopped, statR1unknown, statR2unknown

        recordI2[1] = revcomp(recordI2[1])   #replace with reverse compliment
        if ("N" in recordI1[1] or "N" in recordI2[1]):   #UNKNOWN
            recordI1[1], recordI2[1] = repair(recordI1[1], recordI2[1]) #try to repair N's
            if ("N" in recordI1[1] or "N" in recordI2[1]):  #still have Ns?
                recordR1[0] = recordR1[0] + ": " + recordI1[1] + "-" + recordI2[1] #append index pair to header
                recordR2[0] = recordR2[0] + ": " + recordI1[1] + "-" + recordI2[1]
                R1_unknown.write(f'{recordR1[0]}\n{recordR1[1]}\n{recordR1[2]}\n{recordR1[3]}\n')   #write to file   
                statR1unknown += 1
                R2_unknown.write(f'{recordR2[0]}\n{recordR2[1]}\n{recordR2[2]}\n{recordR2[3]}\n')
                statR2unknown += 1
        if ((recordI1[1] in indexes) and (recordI2[1] in indexes)):     #VALIDATE
            if recordI1[1] == recordI2[1]:                            #MATCHED                                    
                indexpair = recordI1[1] + "-" + recordI2[1]
                indexpermcounts[indexpair] += 1                                          
                recordR1[0] = recordR1[0] + ": " + recordI1[1] + "-" + recordI2[1]
                recordR2[0] = recordR2[0] + ": " + recordI1[1] + "-" + recordI2[1]
                R1_matchfiles[recordI1[1]].write(f'{recordR1[0]}\n{recordR1[1]}\n{recordR1[2]}\n{recordR1[3]}\n')
                statR1match += 1
                R2_matchfiles[recordI2[1]].write(f'{recordR2[0]}\n{recordR2[1]}\n{recordR2[2]}\n{recordR2[3]}\n')
                statR2match += 1
            else: 
                indexpair = recordI1[1] + "-" + recordI2[1]
                indexpermcounts[indexpair] += 1                                                                      #HOPPED
                recordR1[0] = recordR1[0] + ": " + recordI1[1] + "-" + recordI2[1]
                recordR2[0] = recordR2[0] + ": " + recordI1[1] + "-" + recordI2[1]
                R1_hopped.write(f'{recordR1[0]}\n{recordR1[1]}\n{recordR1[2]}\n{recordR1[3]}\n') 
                statR1hopped += 1        
                R2_hopped.write(f'{recordR2[0]}\n{recordR2[1]}\n{recordR2[2]}\n{recordR2[3]}\n')   
                statR2hopped += 1 
        else:
            recordR1[0] = recordR1[0] + ": " + recordI1[1] + "-" + recordI2[1]             #UNKNOWN
            recordR2[0] = recordR2[0] + ": " + recordI1[1] + "-" + recordI2[1]
            R1_unknown.write(f'{recordR1[0]}\n{recordR1[1]}\n{recordR1[2]}\n{recordR1[3]}\n')
            statR1unknown += 1      
            R2_unknown.write(f'{recordR2[0]}\n{recordR2[1]}\n{recordR2[2]}\n{recordR2[3]}\n')  
            statR2unknown += 1

def getstats(indexes: list, statR1match, statR2match, statR1hopped, statR2hopped, statR1unknown, statR2unknown: int) -> dict: #generates plots based on our data
    totallines = statR1match + statR2match + statR1hopped + statR2hopped + statR1unknown + statR2unknown
    totalrecords = totallines/4
    with open("stats.tsv", "a") as statfile:
        statfile.write(f"R1 Matched:\t{statR1match}\n")    #write records to stat file
        statfile.write(f"R2 Matched:\t{statR2match}\n")
        statfile.write(f"R1 Hopped:\t{statR1hopped}\n")
        statfile.write(f"R2 Hopped:\t{statR2hopped}\n")
        statfile.write(f"R1 Unknown:\t{statR1unknown}\n")
        statfile.write(f"R2 Unknown:\t{statR2unknown}\n")
        statfile.write(f"Total Lines:\t{totallines}\n")
        statfile.write(f"Total Records:\t{totalrecords}\n")
        statfile.write(f"\n=Permutations=\n")
        for k,v in indexpermcounts.items():
            statfile.write(f"{k}\t{v}\n")  
    statfile.close()
    records = [statR1match, statR2match, statR1hopped, statR2hopped, statR1unknown, statR2unknown]
    print(records)
    return records

def getplots(records: list): #make bar plot of index-pair permutations
    x = [1,2,3,4,5,6]
    y = records
    xlabels = ["R1 Matched", "R2 Matched", "R1 Hopped", "R2 Hopped", "R1 Unknown", "R2 Unknown"]
    plt.ticklabel_format(style='plain') #no scientific notation
    plt.gcf().subplots_adjust(bottom=0.25)
    plt.bar(x, y, align='center')
    plt.xticks(x, xlabels, rotation=60)
    plt.savefig('counts.png')
    plt.figure().clear()               #clear all plot info from memory (for new plotting)
    plt.close()
    plt.cla()
    plt.clf()
    #start pie chart of records
    y = records
    plt.pie(y, labels=xlabels, autopct='%1.1f%%', shadow=True)
    plt.savefig('records.png')
    return

print("Starting demux...")
statR1match, statR2match, statR1hopped, statR2hopped, statR1unknown, statR2unknown = getrecords(indexes)  #demux engine
print("Generating stats...")
records = getstats(indexes, statR1match, statR2match, statR1hopped, statR2hopped, statR1unknown, statR2unknown)   #get record counts
print("Creating plots...")
getplots(records)             #generate plots from counts
print("Completed.")