#!/usr/bin/env python3

import argparse
from rNMPUtils  import *

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Count context frequency of rNMP incorporation')
    parser.add_argument('GENOME', type=argparse.FileType('r'), help='Reference genome FASTA')
    parser.add_argument('BED', type=argparse.FileType('r'), nargs='+', help='BED file of rNMP incorporation')
    parser.add_argument('-f', action='store_true', help='Use fourth column of bed file as rNMP frequency')
    parser.add_argument('-m', action='store_true', help='Count mononucleotide frequency')
    parser.add_argument('-d', action='store_true', help='Count dinucleotide frequency')
    parser.add_argument('--dist', default=[1], type=int, nargs='+', help='distance between dinucleotides')
    parser.add_argument('-t', action='store_true', help='Count trinucleotide frequency')
    parser.add_argument('-o', default='', help='Output basename')
    args = parser.parse_args()

    if not(any([args.m,args.d, args.t])):
        print('Count mononucleotide by default!')
        args.m = True

    # get position of ribose
    libs, ribos = get_ribo_position(args.BED, args.f)
    print('Ribonucleotides extracted!')

    # get ribos
    results = get_ribo(ribos, libs, args.GENOME, args.m, args.d, args.t, args.dist)
    print('Calculation finished')

    # output
    output(results, args.o)
    print('Done!Output to {}!'.format(args.o))


if __name__ == '__main__':
    main()
