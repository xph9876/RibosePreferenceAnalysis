#!/usr/bin/env python3

import seaborn as sns
import argparse
import matplotlib.pyplot as plt
import pandas as pd


# generate dataframe from tsv file
def generate_df(fr, remove_RE):
    header = fr.readline().rstrip('\n').split('\t')
    header[-1] = 'U'
    data = []
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if len(ws) != len(header):
            continue
        words = ws[0].split('-')
        if remove_RE:
            cate = '-'.join(words[:-2])
            name = words[-2]
        else:
            cate = '-'.join(words[:-1])
            name = words[-1]
        for i in range(1, len(ws)):
            data.append([name, cate, header[i], float(ws[i])])
    df = pd.DataFrame(data, columns=['Library','Genotype', 'rNMP', 'Values'])
    return df

def main():
    parser = argparse.ArgumentParser(description='Generate barplot with datapoints from a mono tsv file')
    parser.add_argument('tsv', type=argparse.FileType('r'), help='Input tsv file')
    parser.add_argument('--legend', action='store_true', help='Draw figure legend on the plots')
    parser.add_argument('--remove_RE', action='store_true', help='Remove RE label in input file')
    parser.add_argument('-o', help='Output plot name')
    args = parser.parse_args()

    if not args.o:
        args.o = args.tsv.name.split('.')[0] + '_bar.png'

    # define color palette
    pal = ['#E31A1C', '#1F78B4', '#FFFFB9', '#33A02C']

    # generate define
    df = generate_df(args.tsv, args.remove_RE)

    # draw
    sns.set(font_scale=2.3, style='ticks')
    fig, ax = plt.subplots(figsize=(8,6))
    plt.subplots_adjust(left=0.1, right=1, top=0.98, bottom=0.2)
    sns.barplot(x='Genotype', y='Values', hue='rNMP', data=df, palette=pal,\
            errorbar='sd', errwidth=1.2, capsize=0.12, edgecolor='k',linewidth=1.5,ax=ax)
    sns.swarmplot(x='Genotype', y='Values', hue='rNMP', data=df, dodge=True,\
            palette='dark:black', ax=ax)
    sns.despine()
    plt.ylim((0,1))
    plt.ylabel('')
    plt.xlabel('')

    # tick labels
    xticklabels = []
    for l in ax.get_xticklabels():
        geno = l.get_text()
        count = len(df[df.Genotype == geno])/4
        xticklabels.append(geno.replace(' ', '\n') + f'\nN = {count:.0f}')
    ax.set_xticklabels(xticklabels,fontsize=19)

    # legend
    if args.legend:
        handles, labels = ax.get_legend_handles_labels()
        plt.legend(handles[4:], labels[4:])
    else:
        ax.get_legend().remove()

    plt.savefig(args.o)

if __name__ == '__main__':
    main()
