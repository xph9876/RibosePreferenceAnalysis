#!/usr/bin/env python3

import argparse
import sys

# check whether the chromosome is demanded
def check_chromosome(s, crs, v):
    if v:
        if s in crs:
            return False
        else:
            return True
    else:
        if s in crs:
            return True
        else:
            return False


def main():
    # argument parser
    parser=argparse.ArgumentParser(description='Find specific chromosome(s) from rNMP incorporation or background frequency files')
    parser.add_argument('infile', nargs='+', type=argparse.FileType('r'), default=sys.stdin, help='rNMP incorporation of background frequency files')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output, default = stdout')
    parser.add_argument('-s', default=['chrM'], nargs = '+', help='Strings to capture, default = chrM')
    parser.add_argument('-v', action='store_true', help='Select non-matching chromosomes')
    parser.add_argument('-a', action='store_true', help='Append to original file')
    parser.add_argument('--name', help='Name for the output line, default = input file name')
    args=parser.parse_args()

    # header
    header = args.infile[0].readline()
    args.o.write(header)
    col_num = len(header.split()) - 1

    # search
    for fr in args.infile:
        fr.seek(0)
        # remove first line
        l = fr.readline()
        freq = []
        for l in fr:
            if args.a:
                args.o.write(l)
            cr = l.split('\t')[0]
            if check_chromosome(cr, args.s, args.v):
                try:
                    freq = [freq[i] + float(l.split('\t')[i+1]) for i in range(col_num)]
                except IndexError:
                    freq = list(map(lambda x : float(x), l.split('\t')[1:]))
        if not args.name:
            name = fr.name.split('/')[-1].split('.')[0]
        else:
            name = args.name
        args.o.write( name +'\t'+ '\t'.join([str(i) for i in freq]) + '\n')

    print("Done!")


if __name__ == '__main__':
    main()
