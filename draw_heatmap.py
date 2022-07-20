#!/usr/bin/env python3

import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools as it
from collections import defaultdict
from matplotlib.ticker import FixedLocator


# read data and check it is mono, di or tri nucleotide
def load_data(data):
    df = pd.read_csv(data, sep='\t').set_index('Sample')
    features = list(df.columns)
    # check mono, di, or tri
    nbase = len(features[0])
    assert len(features) == 4 ** nbase, f"The header line of {data.name} is incorrect."+\
        f"Should be {4**nbase} features. Only found {len(features)} features!"
    # generate feature labels
    rNMPs = 'ACGU'
    dNMPs = 'ACGT'
    if nbase == 1:
        labels = list(rNMPs)
    elif nbase == 2:
        rNMP_loc = 1 - features.index('CA')//4
        labels = list(it.product(rNMPs, dNMPs))
        if rNMP_loc == 0:
            labels = [x[0] + x[1] for x in labels]
        else:
            labels = [x[1] + x[0] for x in labels]
    elif nbase == 3:
        rNMP_loc = 2 - features.index('GCA')//16
        labels = list(it.product(rNMPs, dNMPs, dNMPs))
        if rNMP_loc == 0:
            labels = [x[0] + x[1] + x[2] for x in labels]
        elif rNMP_loc == 1:
            labels = [x[1] + x[0] + x[2] for x in labels]
        else:
            labels = [x[2] + x[1] + x[0] for x in labels]
    return df, labels


# determine the size group from df
def determine_group_size(df):
    data = list(df.iloc[0])
    # add 0.1 to deal with precision issue
    return int(len(data)/sum(data) + 0.1)


# add background frequency to label
def add_bg_freq(labels, group_size, fr, chrom):
    # get bg freq from background
    freqs = {x.replace('U', 'T'):None for x in labels}
    features = fr.readline().rstrip('\n').split('\t')[1:]
    for l in fr:
        ws = l.split('\t')
        if ws[0] == chrom:
            ws = ws[1:]
            for i in range(len(ws)):
                freqs[features[i]] = float(ws[i])
            break
    assert list(freqs.values())[0] != None, f'Cannot find all background frequencies of {chrom} in background frequency file {fr.name}!'
    labels = [[x, freqs[x.replace('U','T')]] for x in labels]
    # sum to 1
    s = []
    for i in range(int(len(labels)/group_size)):
        s += [sum([x[1] for x in labels[i*group_size:i*group_size+group_size]])]*group_size
    for j in range(len(labels)):
        labels[j][1] /= (s[j]/100)
    return [f'{x[1]:.2f}% {x[0]}' for x in labels]



# draw heatmaps
def draw(df, labels, output, no_annot, palette):
    # size parameters
    cell_height = {1:1, 2:0.6, 3:0.3}
    cell_width = {1:1, 2:0.6, 3:0.6}
    font_size = 0.28
    font_sizes_in_cell = {1:20, 2:12, 3:12}

    # get information from df
    nbases = len(df.columns[1])
    samples = list(df.index)

    # set figure size
    longest = max([len(x) for x in samples])
    label_width = len(labels[0]) * font_size + 0.2
    sample_height = longest * font_size
    title_height = 0.3
    colorbar_width = 2
    width = len(samples) * cell_width[nbases] + label_width + colorbar_width 
    height = 4 ** nbases * cell_height[nbases] + sample_height + title_height
    fig, ax = plt.subplots(figsize=(width, height))
    plt.subplots_adjust(left=label_width/width, right=1-colorbar_width/width, \
        top=1-title_height/height, bottom=sample_height/height)
    
    # color scale
    group_size = determine_group_size(df)
    if group_size == 4:
        cmax = 0.5
    elif group_size == 16:
        cmax = 0.125
    
    # draw heatmap
    sns.heatmap(df.T, vmin=0, vmax=cmax, ax=ax, annot=(not no_annot), cmap=palette, annot_kws={"size":font_sizes_in_cell[nbases]})

    # title and axis
    ax.set_xticklabels(samples, rotation='vertical', size=font_size*100)
    ax.set_yticklabels(labels, rotation='horizontal', size=font_size*100)
    ax.set_ylabel('')
    ax.set_xlabel('')

    # change top of colorbar to '0.5-1'
    cax = plt.gcf().axes[-1]
    color_labels = cax.get_ymajorticklabels()
    color_labels_texts = [i.get_text() for i in color_labels]
    color_labels_texts[-1] += ' - 1'
    tick_loc = cax.get_yticks().tolist()
    cax.yaxis.set_major_locator(FixedLocator(tick_loc))
    cax.set_yticklabels(color_labels_texts, fontsize=font_size*100)

    # show or save
    fig.savefig(output)
 

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw heatmap for rNMP incorporation mono-, di-, or tri-nucleotide data.')
    parser.add_argument('DATA',type=argparse.FileType('r'), help='Normalized rNMP incorporation frequncy')
    parser.add_argument('-o', default='rNMP_heatmap.png', help='Output figure name, default= rNMP_heatmap.png')
    parser.add_argument('-b', type=argparse.FileType('r'), help='Select background file. If a file is selected, the background percentage is added to labels.')
    parser.add_argument('--background_chrom', default='chrM', help='Chromosome name of background file, default = chrM, use with -b')
    parser.add_argument('--no_annot', action='store_true', help='Hide percentage annotation in each cell')
    parser.add_argument('--palette', default='icefire', help='Define the palette used for the heatmap, icefire or RdBu_r')
    args = parser.parse_args()


    # set color settings for graph
    sns.set(style='white')

    # get data and information
    df, labels = load_data(args.DATA)
    group_size = determine_group_size(df)
    
    # read background frequency
    if args.b:
        labels = add_bg_freq(labels, group_size, args.b, args.background_chrom)

    # draw heatmaps
    draw(df, labels, args.o, args.no_annot, args.palette)
    print(f'Heatmap is saved to {args.o}!')

if __name__ == '__main__':
    main()
