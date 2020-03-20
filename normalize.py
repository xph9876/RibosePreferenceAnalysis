#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Normalization of nuclear dinucleotide, use get_nuc or get_mito first')
    parser.add_argument('raw', type=argparse.FileType('r'), help='Dinucleotide matrices, raw file')
    parser.add_argument('bg', type=argparse.FileType('r'), help='Background frequency')
    parser.add_argument('-o', type=argparse.FileType('w'),default=sys.stdout, help='Output')
    parser.add_argument('--group_len', default=0, type=int, choices=[0,4,16], help='Number of rows of which the sum is 1, [4,16,0], 0 means all, default=0')
    parser.add_argument('--name', help='Name of chromosome in background frequency used for normalize, default = saccer')
    args = parser.parse_args()

    name = args.bg.readline().rstrip('\n').split()[1:]

    # bg[species][di/tri] = freq
    bg = defaultdict(dict)
    for l in args.bg:
        ws = l.rstrip('\n').split()
        for i in range(1,len(ws)):
            bg[ws[0]][name[i-1]] = float(ws[i])

    # header
    l = args.raw.readline()
    args.o.write(l)
    di = l.rstrip('\n').split('\t')

    # normalize
    for l in args.raw:
        ws = l.split('\t')
        if len(ws) != len(di):
            continue
        freq = []
        # check species
        sp = 'saccer'
        for k in bg.keys():
            if ws[0].find(k)!=-1:
                sp = k
        if args.name:
            sp = args.name
        # normalize freq
        try:
            for i in range(1, len(ws)):
                freq.append(float(ws[i])/bg[sp][di[i]])
        except KeyError:
            sys.exit('[ERROR] Cannot find chrom {} in background file'.format(sp))
        # sum to 1
        s = []
        if args.group_len == 0:
            group_len = len(freq)
        else:
            group_len = args.group_len
        for i in range(int(len(freq)/group_len)):
            s += [sum(freq[i*group_len:i*group_len+group_len])]*group_len
        # prevent 0
        for i in range(len(s)):
            if s[i] == 0:
                s[i] = 2**32
        args.o.write(ws[0] + '\t' + '\t'.join([str(freq[i]/s[i]) for i in range(len(freq))]) + '\n')


    print('Done!')


if __name__ == '__main__':
    main()
