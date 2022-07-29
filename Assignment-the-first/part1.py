import bioinfo
import numpy as np
import matplotlib.pyplot as plt
import gzip
import argparse as ap

parser = ap.ArgumentParser(description='Generate mean scores and plot for 101bp FQ reads')
parser.add_argument('-i', '--infile', help='In-file (fq.gz)', type=str)
parser.add_argument('-o', '--outplot', help='Out-Plot (.png)', type=str)
parser.add_argument('-s', '--seqlen', help='Sequence Length', type=int) 
args = parser.parse_args()
inFile = args.infile
outPlot = args.outplot
seqLen = args.seqlen


def getscores(inFile):
    scores = np.zeros(int(seqLen))
    numlines = 0
    with gzip.open(inFile, 'rt', encoding='utf-8') as fileBio1:
        for line in fileBio1:
            numlines += 1
            if numlines % 4 == 0:
                line=line.strip()
                for k, v in enumerate(line):
                    scores[k] += bioinfo.convert_phred(v)
        for k, v in enumerate(scores):
            scores[k] = scores[k]/(numlines/4)
    return(scores)

scores = getscores(inFile)
print(scores)
x = list(range(1,(seqLen+1)))
y = scores
plt.plot(x,y)
plt.xlabel('Position')
plt.ylabel('Mean Score')
plt.title('Mean scores over sequence position')
plt.savefig(outPlot)