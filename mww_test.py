#!/usr/bin/env python3

import sys
import argparse
from scipy.stats import mannwhitneyu as mww


# read data
def read_data(fr):
    data = {}
    header = fr.readline().rstrip('\n').split('\t')[1:]
    data['all'] = {k:[] for k in header}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        name = ws[0]
        c = name.split('-')
        cate = '-'.join(c[:-1])
        if cate not in data:
            data[cate] = {k:[] for k in header}
        for i in range(len(ws[1:])):
            data[cate][header[i]].append(float(ws[i+1]))
            data['all'][header[i]].append(float(ws[i+1]))
    return data, header


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Perform Mann-Whitney U test for heatmap data')
    parser.add_argument('tsv', type=argparse.FileType('r'), help='Normalized frequency file')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('--baseline', type=float, default=-1, help='Baseline value to compare with, default: infer from heatmap data')
    parser.add_argument('--min', type=int, default=4, help='Min sample number to keep the group, default: 4')
    args = parser.parse_args()

    # read data
    data, header = read_data(args.tsv)
    n = len(header[0])
    if args.baseline == -1:
        if n == 1 or n == 2:
            args.baseline = 0.25
        else:
            args.baseline = 1/16

    # perform mww test
    args.o.write(f'Celltype\t' + '\t'.join(header) + '\n')
    for celltype, v in data.items():
        # check
        if len(list(v.values())[0]) < args.min:
            continue
        args.o.write(f'{celltype}')
        for k in header:
            d = v[k]
            _, p = mww(d, [args.baseline]*len(d), alternative='greater')
            args.o.write(f'\t{p}')
        args.o.write('\n')
        
    print('Done!')

if __name__ == '__main__':
    main()
