# Compare ID in bnx file to ID in Xmap file
# Updated - Checks for molecules that are 150kb or larger &
# have sites within 2.5kb of each end that have SNR > 50 & Label Intensity > 5
#Output is a bnx file with only telomeres and a text file with list of IDs
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

def getmatchedSNRsites(SNRs, min_snr):
    SNR_site = 0
    matchedSNRs = []
    for SNR in SNRs:
        SNR_site = SNR_site + 1
        if float(SNR) >= min_snr:
            matchedSNRs.append(SNR_site)
    return matchedSNRs


def getmatchedIntensities(intensities, matchedSNRs, min_intensity):
    inten_site = 0
    matchedInten = []
    for inten in intensities:
        inten_site = inten_site + 1
        if inten_site in matchedSNRs and float(inten) >= min_intensity:
            matchedInten.append(inten_site)
    return matchedInten

def main():
    min_mol_length = 150000
    min_distance = 3000
    min_snr = 50
    min_intensity = 3.5
    min_sites = 4

    bnxfile = 'all.bnx'
    filename = bnxfile.split('/')[-1].split('.')[0]
    outbnx = filename + '_telomeres.bnx'
    outidlist = filename + '_telomere_ids.txt'

    matched = []
    lines = []
    header = []
    with open(bnxfile, 'r') as f:
        while True:
            line = f.readline()
            if line.startswith('#'):
                header.append(line)
            if not line.startswith('#'):
                lines.append(line)
                break
        lines = lines + f.readlines()


    match = False
    for line in lines:
        #checks length of molecule
        if line.startswith('0'):
            match = False
            id = str(line.split('\t')[1])
            length = float(line.split('\t')[2])
            # print(length)
            if length >= min_mol_length:
                match = True
                tel = []
                # print(line)
                tel.append(line)
            else:
                match = False
        #picks sites within 2500kb
        elif match and line.startswith('1'):
            left_count = 0
            right_count = 0
            sites = line.split('\t')[1:]
            if len(sites) >= min_sites:
                for site in sites:
                    if float(site) < min_distance:
                        # print(site)
                        left_count = left_count + 1
                    elif float(site) > (length-min_distance):
                        # print(site)
                        right_count = right_count + 1
                right_count = right_count - 1
                if left_count > 0 or right_count > 0:
                    tel.append(line)
                else:
                    tel = []
                    match = False
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
                leftmatchedSNRs = getmatchedSNRsites(SNRsleft, min_snr)
                # print(leftmatchedSNRs)
            #if possible telomere on right
            if right_count > 0:
                right_count = right_count * -1
                SNRsright = line.split('\t')[right_count:]
                rightmatchedSNRs = getmatchedSNRsites(SNRsright, min_snr)
                # print(rightmatchedSNRs)
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
                leftmatchedInten = getmatchedIntensities(intensities, leftmatchedSNRs, min_intensity)
            # if possible telomere on right
            if len(leftmatchedInten) < 1 and len(rightmatchedSNRs) > 0:
                intensities = line.split('\t')[right_count:]
                rightmatchedInten = getmatchedIntensities(intensities, rightmatchedSNRs, min_intensity)
                # print(rightmatchedInten)
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

    bnxids = []
    with open(outbnx, 'w') as o:
        for h in header:
            o.write(h)
        for line in flat_list:
            o.write(line)


    with open(outbnx, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        idlines = [line for line in reader if '#' not in line[0]]
        bnxids = [row[1] for row in idlines if '0' in row[0]]
        while "" in bnxids:
            bnxids.remove("")

    with open(outidlist, 'w') as o:
        for id in bnxids:
            o.write(id + '\n')


if __name__ == '__main__':
    main()