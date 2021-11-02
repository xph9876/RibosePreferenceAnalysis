#!/usr/bin/env python3

import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import itertools as it
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


# draw color bar
def draw(df, labels, output, palette, font_size):
    fig, ax = plt.subplots(figsize=(6, 1))
    plt.subplots_adjust(bottom=0.5)
    
    # color scale
    group_size = determine_group_size(df)
    if group_size == 4:
        cmax = 0.5
    elif group_size == 16:
        cmax = 0.125
    
    # draw heatmap
    cmap = sns.color_palette(palette, as_cmap=True)
    norm = mpl.colors.Normalize(vmin=0, vmax=cmax)
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal')

    # # change top of colorbar to '0.5-1'
    ax.xaxis.set_ticks(np.arange(0, cmax*1.01, cmax/5))
    fig.canvas.draw()
    labels = ax.get_xmajorticklabels()
    labels_texts = [i.get_text() for i in labels]
    labels_texts[-1] += ' - 1'
    tick_loc = ax.get_xticks().tolist()
    ax.xaxis.set_major_locator(FixedLocator(tick_loc))
    ax.set_xticklabels(labels_texts, fontsize=font_size)

    # show or save
    plt.savefig(output)
 

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw color scale for the heatmaps.')
    parser.add_argument('DATA',type=argparse.FileType('r'), help='Normalized rNMP incorporation frequncy')
    parser.add_argument('-o', default='cbar.png', help='Output figure name, default= rNMP_heatmap.png')
    parser.add_argument('--palette', default='icefire', help='Seaborn palette used for the heatmap')
    parser.add_argument('--font_size', default=18, type=float, help='Font scale for the color bar')
    args = parser.parse_args()

    # get data and information
    df, labels = load_data(args.DATA)
    group_size = determine_group_size(df)
    
    # draw heatmaps
    draw(df, labels, args.o, args.palette, args.font_size)
    print(f'Color scale plot is saved to {args.o}!')

if __name__ == '__main__':
    main()