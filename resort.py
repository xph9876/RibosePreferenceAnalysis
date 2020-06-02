#!/usr/bin/env python3

import argparse
import sys

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Resort to specific order')
    parser.add_argument('infile', type=argparse.FileType('r'), help='File needed to resort')
    parser.add_argument('order', type=argparse.FileType('r'), help='Order needed')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output file')
    parser.add_argument('-d', default='-', help='Connector of library information, default = \'-\'')
    parser.add_argument('-c', type=int, default=1, help='Column number of library num, default=1')

    args = parser.parse_args()
    args.c -= 1

    # get header
    args.o.write(args.infile.readline())

    # resort
    data = {}
    for l in args.infile:
        fs,d = l.split('\t',1)
        data[fs] = d

    # get fs order
    for l in args.order:
        ws = l.rstrip('\n').split('\t')
        fs = ws[args.c]
        if fs in data:
            args.o.write('{}\t{}'.format(args.d.join(ws), data[fs]))

    print('Done!')


if __name__ == '__main__':
    main()
