#!/usr/bin/env python3

import argparse
import sys

def rc(s):
    r = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    out = ''
    for c in s:
        out = r[c] + out
    return out

def main():
    parser = argparse.ArgumentParser(description='Generate background file for the opposite strand')
    parser.add_argument('tsv', type=argparse.FileType('r'), help='TSV file of background frequency')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file')
    args = parser.parse_args()

    # load data
    header = args.tsv.readline()
    features = header.rstrip('\n').split('\t')
    args.o.write(header)
    for l in args.tsv:
        ws = l.rstrip('\n').split('\t')
        data = {rc(features[i]):ws[i] for i in range(1, len(features))}
        if ws[0].endswith(r'-)'):
            ws[0] = ws[0].replace(r'(-)', r'(+)')
        elif ws[0].endswith(r'+)'):
            ws[0] = ws[0].replace(r'(+)', r'(-)')
        data[features[0]] = ws[0]
        args.o.write('\t'.join([data[x] for x in features]) + '\n')
    
    print('Done')


if __name__ == '__main__':
    main()
