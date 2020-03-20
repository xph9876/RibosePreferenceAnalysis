#!/usr/bin/env python3

from collections import defaultdict
import argparse
import sys

# reversecomplement
def rc(s):
    c = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    r = ''
    for i in s:
        r = c[i] + r
    return r

def main():
    # argument parser
    parser = argparse.ArgumentParser(description='Count dinucleotides from fasta file')
    parser.add_argument('FASTA', type=argparse.FileType('r'), default=sys.stdin, help='Input fasta file to count')
    parser.add_argument('-d', type=int, default=1, help='Distance between dinucleotide pair, default=1')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('-s', action='store_true', help='Count only one strand')
    parser.add_argument('--mono', action='store_true', help='Count mono nucleotide instead')
    parser.add_argument('--trinuc', action='store_true', help='Count trinucleotide instead')
    args = parser.parse_args()

    if args.mono and args.trinuc:
        sys.exit('[ERROR] Can only count mononucleotide or trinucleotide!')

    # count
    data = defaultdict(lambda : defaultdict(int))
    cache = ''
    for l in args.FASTA:
        l = l.rstrip('\n')
        if l[0] == '>':
            chrom = l[1:]
            print('Start counting {}!'.format(chrom))
            cache = ''
        elif args.mono:
            for i in l:
                data[chrom][i] += 1
        elif args.trinuc:
            l = cache + l
            for i in range(len(l)-2):
                data[chrom][l[i:i+3]] += 1
            cache = l[-2]
        else:
            l = cache + l
            for i in range(len(l)-args.d):
                data[chrom][l[i] + l[i+args.d]] += 1
            cache = l[-(args.d):]

    # summarize
    result = defaultdict(lambda : defaultdict(int))
    for chrom, v in data.items():
        for kmer, count in v.items():
            result[chrom][kmer.upper()] += count

    # generate nuc
    base = ['A','C','G','T']
    if args.mono:
        nuc = base
    elif args.trinuc:
        nuc = []
        for i in base:
            for j in base:
                for k in base:
                    nuc.append(i+j+k)
    else:
        nuc = []
        for i in base:
            for j in base:
                nuc.append(i+j)

    # single strand
    if not args.s:
        for chrom, v in data.items():
            for n, count in v.items():
                nu = n.upper()
                if nu in nuc:
                    result[chrom][rc(nu)] += count

    # output
    args.o.write('chr\t' + '\t'.join(nuc) + '\n')
    for k in sorted(result.keys()):
        args.o.write('{}\t'.format(k) + '\t'.join([str(result[k][i]) for i in nuc]) + '\n')
    print('Done!')

if __name__ == '__main__':
    main()

