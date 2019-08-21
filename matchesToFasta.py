import pandas as pd
import sys


def search(arr, s):
    counter = 0
    for j in range(len(arr)):

        # checking if string given in query
        # is present in the given string.
        # If present, increase times
        if (s == (arr[j])):
            counter += 1

    return counter


def main(matchfile, outfasta):
# def main():
    # matchfile = './genome/allchroms.txt'
    sko = pd.read_csv(matchfile, sep='\t', header=0, index_col=None)
    #print(sko.head(5))
    outputchr = matchfile.split('.')[0]
    # outputpng = outputchr + '_NEWhist.png'
    # outputdf = outputchr + '_NEWfrequency.tsv'
    reads = []
    smallseq = ''
    matchdict = {}
    for index, row in sko.iterrows():
        read = row['seqID']
        #print(read)
        if read not in reads:
            reads.append(read)
            smallseq = ''
            match = str(row['matched']).upper()[-1]
            smallseq = smallseq + match
            matchdict[read] = smallseq
        else:
            match = str(row['matched']).upper()[-1]
            smallseq = smallseq + match
            matchdict[read] = smallseq
        #print(len(reads))
    print('small seq found')
    matches = pd.DataFrame.from_dict(matchdict, orient='index')
    # print(matches)
    # matches.to_csv('allchroms.tsv', sep='\t', header=None)
    with open(outfasta, 'a') as o:
        for index, row in matches.iterrows():
            name = index
            seq = row[0]
            string = '>' + str(name) + '\n' + seq + '\n'
            o.write(string)

if __name__ == '__main__':
    matchfile = sys.argv[1]
    outfasta = sys.argv[2]
    main(matchfile, outfasta)
    # main()
