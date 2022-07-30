# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Bio1  | 101 | Phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | Index1 | 8 | Phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | Index2 | 8 | Phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | Bio2 | 101 | Phred+33 |
```
zcat (file) | head -8   # for each file
zcat (file) | head -2 | tail -1 | wc 
```


2. Per-base NT distribution


![R1](plotR1.png)
![R2](plotR2.png)
![R3](plotR3.png)
![R4](plotR4.png)


   3. For R1, the reads were good so 36 seems like a good cutoff. If anything fell below 36 for no good reason, I would not include it. For example, the read quality diminishes over time and gets close to 36 towards the end of the read (position 90+) and the last couple of nucleotides drop below 36 in R4. For index 1 and index 2, they increase over time and are all over 30, so they aren't terrible, but I would say 36 is a goal to obtain. However, given the data we have, we could potentially use all of the index data by finding confidence scores between an index pair. For example: NAGTTCATC and CAGTTCATN, we could fill in the N-ends with its sister index's end and be pretty confident that they are the same index. 
    
   4. zcat "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" | sed -n '/^@/{n;p;}' | awk '/N/ {count++} END{print count}'
R2: 3976613
R3: 3328051
Total: 7304664
    
## Part 2
1. Define the problem
```
Reads are contained in the file ordered by X,Y location from the flowcell and not by indexes.
All indexes are included, even if they are bad reads or hopped.
All reads, regardless of index, are in 2 files (R1 and R4)
```
2. Describe output
```
We want to output:
   the reads with their index pair in the header mapped to an seperate read file for each index (24 * 2)
   the hopped/unmatched indexes in to a seperate file (2 files)
   the unknown indexes in to a seperate file (2 files)
```
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
```
Done
```
4. Pseudocode
```
Done
```
5. High level functions. For each function, be sure to include:
    1. Description/doc string
```
Done
```
   2. Function headers (name and parameters)
```
Done
```
   3. Test examples for individual functions
```
Done
```
   4. Return statement
```
Done
```
