#select 50 random ranges of 1Mb and write them to .fa files

import numpy as np
import sys

def main(chrfile):
    # chrfile = 'chr22.fa'
    chr = chrfile.split('/')[-1].split('.')[0]
    outputfile = chr + 'reads.fa'

    with open(chrfile, 'r') as f:
        lines = f.read().splitlines()

    sequence = ''.join(lines)
    lastbase = len(sequence)
    readnum = round((30000*lastbase)/3100000000)

    windows = []
    for i in range(readnum):
        start = np.random.randint(lastbase)
        end = start+1000001
        tup = (start, end)
        windows.append(tup)

    with open(outputfile, 'a') as o:
        for window in windows:
            st = window[0]
            e = window[1]
            seq = sequence[st:e]
            string = '>' + str(st) + '-' + str(e) + '\n' + seq + '\n'
            o.write(string)





if __name__ == '__main__':
    chrfile = sys.argv[1]
    main(chrfile)