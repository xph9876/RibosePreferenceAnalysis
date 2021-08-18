#!/usr/bin/env python3

import argparse
import sys

# load bg
def load_bg(fr):
    fea = fr.readline().rstrip('\n').split('\t')
    bg = {}
    for l in fr:
        ws =l.rstrip('\n').split('\t')
        if len(ws) > 1:
            bg[ws[0]] = {}
            for i in range(1, len(ws)):
                bg[ws[0]][fea[i]] = float(ws[i])
    return bg


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Normalization of rNMP preference analysis')
    parser.add_argument('raw', type=argparse.FileType('r'), help='rNMP incorporation file for desired chromosome.')
    parser.add_argument('bg', type=argparse.FileType('r'), help='Background frequency')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file.')
    parser.add_argument('--group_len', default=0, type=int, choices=[0,4,16], help='Number of rows of which the sum is 1, [4,16,0]." + \
                        " if 0 is selected, the sum of all rows would be 1. default = 0')
    args = parser.parse_args()

    # load bg
    bg = load_bg(args.bg)

    # header
    l = args.raw.readline()
    args.o.write(l)
    fea = l.rstrip('\n').split('\t')

    # normalize
    for l in args.raw:
        ws = l.split('\t')
        if len(ws) != len(fea):
            continue
        # normalize freq
        freq = []
        sp = ws[0]
        try:
            for i in range(1, len(ws)):
                freq.append(float(ws[i])/bg[sp][fea[i]])
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
    print ('Done!')


if __name__ == '__main__':
    main()