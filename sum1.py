#!/usr/bin/env python3

import argparse
import sys

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Change the sum of each ribos to 1')
    parser.add_argument('infile', type=argparse.FileType('r'), nargs='?', default=sys.stdin, help='infile')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output file')
    args = parser.parse_args()

    # read
    l = args.infile.readline()
    args.o.write(l)
    for l in args.infile:
        name = l.split('\t')[0]
        freq = [float(i) for i in l.rstrip('\n').split('\t')[1:]]
        s = []
        for i in range(int(len(freq)/4)):
            s += [sum(freq[i*4:i*4+4])]*4
        args.o.write('\t'.join([name]+[str(freq[i]/s[i]) for i in range(len(freq))]) +'\n')

    print('Done!')

if __name__ == '__main__':
    main()
