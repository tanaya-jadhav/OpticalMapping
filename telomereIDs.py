# Compare ID in bnx file to ID in Xmap file
# Updated - Checks for molecules that are 150kb or larger &
# have sites within 2.5kb of each end that have SNR > 50 & Label Intensity > 5
# author: Tanaya Jadhav

import csv


def idlist(file, column):
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        lines = [line for line in reader if '#' not in line[0]]
        ids = [row[column] for row in lines]
        while "" in ids:
            ids.remove("")
        return ids


def getmatchedSNRsites(SNRs):
    SNR_site = 0
    matchedSNRs = []
    for SNR in SNRs:
        SNR_site = SNR_site + 1
        if float(SNR) >= 50:
            matchedSNRs.append(SNR_site)
    return matchedSNRs


def getmatchedIntensities(intensities, matchedSNRs):
    inten_site = 0
    matchedInten = []
    for inten in intensities:
        inten_site = inten_site + 1
        if inten_site in matchedSNRs and float(inten) > 5:
            matchedInten.append(inten_site)
    return matchedInten

def main():
    bnxfile = 'SKMEL2_filtered_merged.bnx'
    # xmapfile = 'EXP_REFINEFINAL1_chromosome16_Filtered.xmap'
    # xmapIDs = idlist(xmapfile, 1)

    matched = []
    lines = []
    with open(bnxfile, 'r') as f:
        while True:
            line = f.readline()
            if not line.startswith('#'):
                lines.append(line)
                break

        lines = lines + f.readlines()
    match = False
    for line in lines:
        #checks length of molecule
        if line.startswith('0'):
            match = False
            length = float(line.split('\t')[2])
            # print(length)
            if length > 150000:
                match = True
                tel = []
                tel.append(line)
            else:
                match = False
        #picks sites within 2500kb
        elif match and line.startswith('1'):
            left_count = 0
            right_count = 0
            sites = line.split('\t')[1:]
            for site in sites:
                if float(site) < 2500:
                    left_count = left_count + 1
                elif float(site) > (length-2500):
                    right_count = right_count + 1
            if left_count > 0 or right_count > 0:
                tel.append(line)
            else:
                tel = []
                match = False

        #picks sites with Label SNR >= 50
        elif match and line.startswith('QX11'):
            leftmatchedSNRs = []
            rightmatchedSNRs = []
            #if possible telomere on left
            if left_count > 0:
                SNRsleft = line.split('\t')[1:left_count+1]
                leftmatchedSNRs = getmatchedSNRsites(SNRsleft)
            #if possible telomere on right
            elif right_count > 0:
                right_count = right_count * -1
                SNRsright = line.split('\t')[right_count:]
                rightmatchedSNRs = getmatchedSNRsites(SNRsright)
            if len(leftmatchedSNRs) > 0 or len(rightmatchedSNRs) > 0:
                tel.append(line)
            else:
                tel = []
                match = False

        #picks sites with Label Intensity > 5
        elif match and line.startswith('QX12'):
            leftmatchedInten = []
            rightmatchedInten = []
            # if possible telomere on left
            if len(leftmatchedSNRs) > 0:
                intensities = line.split('\t')[1:left_count+1]
                leftmatchedInten = getmatchedIntensities(intensities, leftmatchedSNRs)
            # if possible telomere on right
            elif len(rightmatchedSNRs) > 0:
                intensities = line.split('\t')[right_count:]
                rightmatchedInten = getmatchedIntensities(intensities, rightmatchedSNRs)

            if len(leftmatchedInten) > 0 or len(rightmatchedInten) > 0:
                tel.append(line)
                matched.append(tel)
            else:
                tel = []
                match = False
        #write rest of the lines of a picked molecule
        elif match:
            matched.append(line)

    flat_list = [item for sublist in matched for item in sublist]

    with open('outputtestboth_new.txt', 'w') as o:
        for line in flat_list:
            o.write(line)


if __name__ == '__main__':
    main()