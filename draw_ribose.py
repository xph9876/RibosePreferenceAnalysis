#!/usr/bin/env python3

import numpy as np
import argparse
import sys
# from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from collections import defaultdict

# argparse
parser = argparse.ArgumentParser(description='PCA for dinucleotide data')
parser_io = parser.add_argument_group('I/O')
parser_io.add_argument('DATA',type=argparse.FileType('r'), default=sys.stdin, help='Normalized dinucleotide data')
parser_io.add_argument('-o', default='', help='output basename')
parser_io.add_argument('-b', type=argparse.FileType('r'), help='Select background file. If a file is selected, the background percentage is added to labels.')
parser_io.add_argument('--nr', action='store_true', help='Input file is NR type dinucleotide frequency. Of which the second position is rNMP.')
parser_io.add_argument('--mono', action='store_true', help='Input file is mononucleotide frequency.')
parser_io.add_argument('--tri', type=int, default=0, help='Input file is trinucleotide frequency, with rNMP incorporated at TRI position')
parser_io.add_argument('--background_chrom', default='chrM', help='Chromosome name of background file, default = chrM')
parser_graphs = parser.add_argument_group('Graphs')
parser_graphs.add_argument('-m', action='store_true', help='Draw heatmap')
parser_graphs.add_argument('-d', action='store_true', help='Draw dendrogram for cluster')
parser_graphs.add_argument('-p', action='store_true', help='Draw PCA scatter plot')
parser_graphs.add_argument('-s', action='store_true', help='Draw Singular value plot')
parser_h = parser.add_argument_group('Heatmap')
parser_h.add_argument('--legend_group', default=0, choices=[0,4,16], type=int, help='Number of lables of which the sum is 1. If 0 is selected, the sum of all labels will be 1. default = 0.')
parser_h.add_argument('--no_annot', action='store_true', help='Hide percentage annotation in each cell')
parser_h.add_argument('--cmax', type=float, default=0.5, help='Maximum value in color scale. Any preferency beyond that will show as the maximum color.')
parser_d = parser.add_argument_group('Cluster dendrogram')
parser_d.add_argument('--chrom_num', dest='cn', type=int, default=18, help='Number of chromosomes, default=18')
parser_p = parser.add_argument_group('PCA graph')
parser_p.add_argument('--3d', dest='td', action='store_true', help='Generate 3D graph')
args = parser.parse_args()

# argument relations
if any([args.m, args.d, args.p, args.s]) == False:
    args.m = True

if sum([args.nr, args.mono, args.tri!=0]) > 1:
    print('Arguments conflict, can only be mono, di or trinuc')
    sys.exit(1)

if args.mono:
    args.legend_sum1 = True


# build ndarray
sample_raw = []
da = []
di = args.DATA.readline().rstrip('\n').split('\t')[1:]
for l in args.DATA:
    ws = l.rstrip('\n').split('\t')
    if len(ws) < len(di):
        break
    try:
        da.append([float(i) for i in ws[1:]])
    except ValueError:
        print(ws)
        sys.exit(1)
    sample_raw.append(ws[0])
mat = np.asarray(da)
sample = [i.split('/')[-1] for i in sample_raw]

# set color settings for graph
sns.set(palette='muted', style='white')

# heatmap
# get percentage from background
labels=[]
f = defaultdict(float)
if args.b:
    for l in args.b:
        ws = l.split('\t')
        if ws[0] == args.background_chrom:
            ws = ws[1:]
            for i in range(len(ws)):
                f[di[i]] += float(ws[i])
else:
    for d in di:
        f[d] = 0
for k,v in f.items():
    labels.append([k,v])

# sort labels
if args.mono:
    labels = sorted(labels, key=lambda x:x[0])
    labels[-1][0]='U'

elif args.tri:
    base_order = [[0,1,2], [1,0,2],[2,1,0]]
    labels = sorted(labels, key=lambda x:[x[0][i] for i in base_order[args.tri-1]])
    # change T to U
    for i in range(len(labels)-16, len(labels)):
        labels[i][0] = labels[i][0][0:args.tri-1] + 'U' + labels[i][0][args.tri:]

elif args.nr:
    labels = sorted(labels, key=lambda x:(x[0][1],x[0][0]))
    labels[-4][0]='AU'
    labels[-3][0]='CU'
    labels[-2][0]='GU'
    labels[-1][0]='TU'
else:
    labels = sorted(labels, key=lambda x:(x[0][0],x[0][1]))
    labels[-4][0]='UA'
    labels[-3][0]='UC'
    labels[-2][0]='UG'
    labels[-1][0]='UT'

if args.b:
    # sum to 1
    if args.legend_group == 0:
        group_len = len(labels)
    else:
        group_len = args.legend_group
    s = []
    for i in range(int(len(labels)/group_len)):
        s += [sum([x[1] for x in labels[i*group_len:i*group_len+group_len]])]*group_len
    for j in range(len(labels)):
        labels[j][1] /= (s[j]/100)

# merge labels
label_texts = []
for i in labels:
    if args.b:
        label_texts.append('{:.2f}% {}'.format(i[1], i[0]))
    else:
        label_texts.append(i[0])


if args.m:
    if args.mono:
        fig, ax = plt.subplots(figsize=(mat.shape[0]*0.45+6,6))
        plt.subplots_adjust(left=0.03, right=1, bottom=0.4)
        sns.heatmap(mat.T, vmin=0, vmax=args.cmax,center=args.cmax*0.45, ax=ax, annot=(not args.no_annot), annot_kws={"size":12})
    elif args.tri:
        fig, ax = plt.subplots(figsize=(mat.shape[0]*0.45+8,16))
        plt.subplots_adjust(left=0.03, right=1, bottom=0.15, top=0.98)
        sns.heatmap(mat.T, vmin=0, vmax=args.cmax,center=args.cmax*0.45, ax=ax, yticklabels=1, annot=(not args.no_annot), annot_kws={"size":12})
    else:
        fig, ax = plt.subplots(figsize=(mat.shape[0]*0.45+6,12))
        plt.subplots_adjust(left=0.03, right=1, bottom=0.3, top=0.95)
        sns.heatmap(mat.T, vmin=0, vmax=args.cmax,center=args.cmax*0.45, ax=ax, annot=(not args.no_annot), annot_kws={"size":12})
    ax.set_title('Heatmap for dinucleotides')
    ax.set_xticklabels(sample, rotation='vertical')
    ax.set_yticklabels(label_texts,rotation='horizontal')

    # change top of colorbar to '0.5-1'
    cax = plt.gcf().axes[-1]
    color_labels = cax.get_ymajorticklabels()
    color_labels_texts = [i.get_text() for i in color_labels]
    color_labels_texts[-1] += ' - 1'
    cax.set_yticklabels(color_labels_texts)

    # add percentage labels
    if args.b:
        ax.set_yticklabels(label_texts)
    # show or save
    if args.o == '':
        plt.show()
    else:
        ax.get_figure().savefig(args.o + '_heatmap.png')
        print('Heatmap is saved to {}_heatmap.png'.format(args.o))

# hierarchy cluster
if args.d:
    fig, ax = plt.subplots(figsize=(mat.shape[0]*0.2+2,8))
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.3)
    z = hierarchy.linkage(mat,method='average')
    dn = hierarchy.dendrogram(z, labels=sample, ax=ax)
    ax.set_title('Cluster for dinucleotides')
    ax.set_ylabel('Distance')

    # set color
    colors = ['darkred','red', 'darkblue', 'blue', 'darkgreen', 'green',\
            'black', 'grey', 'darkpurple', 'purple', '#d7df00', '#fdbf6f']
    color_dic = {}
    for i in range(len(sample)):
        color_dic[sample[i]] = colors[int(i/args.cn)]
    labels = ax.get_xmajorticklabels()
    for i in labels:
        # i.set_color(color_dic[i.get_text()])
        i.set_fontsize(10)
        i.set_rotation('vertical')

    # show or save
    if args.o == '':
        plt.show()
    else:
        ax.get_figure().savefig(args.o + '_dendrogram.png')
        print('Cluster dendrogram is saved to {}_dendrogram.png'.format(args.o))


# Singular value
if args.s:
    # sigular value
    u,s,v = np.linalg.svd(mat)
    fig, ax = plt.subplots()
    plt.plot(range(1,len(s)+1), s, marker='o', ax=ax)
    ax.set_title('Singular value of dinucleotides')
    ax.set_xlabel('Number')
    ax.set_ylabel('Value')

    # show or save
    if args.o == '':
        plt.show()
    else:
        ax.get_figure().savefig(args.o + '_sv.png')
        print('Singular value graph is saved to {}_sv.png'.format(args.o))

# PCA
# if args.p:
    # if args.td:
        # pca = PCA(n_components=3)
    # else:
        # pca = PCA(n_components=2)
    # projected = pca.fit_transform(mat)
    # colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']
    # if args.td:
        # ax = plt.subplot(111, projection='3d')
    # else:
        # fig,ax = plt.subplots(figsize=(10,8))


    # # add dots
    # dots = []
    # start = 0
    # for s in range(len(groups)):
        # data_ziped = list(zip(*projected[start:start + group_num[s]]))
        # if args.td:
            # dots.append(ax.scatter(data_ziped[0], data_ziped[1], data_ziped[2], c=colors[s], s=30))
        # else:
            # dots.append(ax.scatter(data_ziped[0], data_ziped[1], c=colors[s], s=30))
            # for i in range(len(data_ziped[0])):
                # ax.text(data_ziped[0][i], data_ziped[1][i], sample[start+i], fontsize=8, color='grey')
        # start += group_num[s]

    # # add legend
    # ax.legend(dots,groups, scatterpoints=1, fontsize=8, ncol=2)
    # ax.set_title('PCA for dinucleotides')
    # ax.set_xlabel('Component 1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0]*100))
    # ax.set_ylabel('Component 2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1]*100))
    # if args.td:
        # ax.set_zlabel('Component 3 ({:.2f}%)'.format(pca.explained_variance_ratio_[2]*100))

    # # show or save
    # if args.o == '':
        # plt.show()
    # else:
        # ax.get_figure().savefig(args.o + '_cluster.png')
#         print('PCA cluster graph is saved to {}_cluster.png'.format(args.o))
