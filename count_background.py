#!/usr/bin/env python3

from collections import defaultdict
import argparse
import sys
import itertools as it

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
    parser.add_argument('--allow_dup_chroms', action='store_true', help='Sum up all counts of chromosomes with the same name. By default, only the first one is counted')
    args = parser.parse_args()

    if args.mono and args.trinuc:
        sys.exit('[ERROR] Can only count mono-, dinuc- or tri-nucleotide!')

    # count
    data = defaultdict(lambda : defaultdict(int))
    cache = ''
    seen = set()
    for l in args.FASTA:
        l = l.rstrip('\n')
        if l[0] == '>':
            chrom = l[1:]
            if chrom in seen and not args.allow_dup_chroms:
                print('Skip duplicate chromosome {}!'.format(chrom))
                chrom = None
            else:
                seen.add(chrom)
                print('Start counting {}!'.format(chrom))
            cache = ''
        elif args.mono:
            for i in l:
                data[chrom][i] += 1
        elif args.trinuc:
            l = cache + l
            for i in range(len(l)-2):
                data[chrom][l[i:i+3]] += 1
            cache = l[-2:]
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

    # generate kmer list
    bases = 'ACGT'
    if args.mono:
        nuc = list(bases)
    elif args.trinuc:
        nuc = [''.join(x) for x in it.product(bases, repeat=3)]
    else:
        nuc = [''.join(x) for x in it.product(bases, repeat=2)]

    # Add opposite strand for both strands
    if not args.s:
        for chrom, v in data.items():
            for n, count in v.items():
                if n.upper() in nuc:
                    result[chrom][rc(n.upper())] += count
    
    # remove skipped choromosomes
    del result[None]

    # output
    args.o.write('Sample\t' + '\t'.join(nuc) + '\n')
    for chrom in sorted(result.keys()):
        args.o.write(f'{chrom}\t' + '\t'.join([str(result[chrom][kmer]) for kmer in nuc]) + '\n')
    print('Done!')

if __name__ == '__main__':
    main()

