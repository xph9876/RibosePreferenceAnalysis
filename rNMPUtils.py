from collections import defaultdict

# complement
def rc(s):
    c = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    r = ''
    for i in s:
        r = c[i] + r
    return r


# get all positions
def get_ribo_position(frs, use_frequency):
    libs = []
    pos = defaultdict(lambda : [0]*len(frs))
    for fridx in range(len(frs)):
        fr = frs[fridx]
        libs.append(fr.name.split('/')[-1].split('.')[0])
        for l in fr:
            ws = l.rstrip('\n').split('\t')
            if len(ws) != 6:
                continue
            if use_frequency:
                pos[(ws[0], int(ws[2]), ws[5])][fridx] += float(ws[3])
            else:
                pos[(ws[0], int(ws[2]), ws[5])][fridx] += 1

    results = defaultdict(list)
    for k in sorted(pos.keys()):
        results[k[0]].append([k[1], k[2], pos[k]])
        del pos[k]
    return libs, results


# read genome and count
def get_ribo(ribos, libs, gr, mono, dinuc, trinuc, dist):
    # initial output
    result = {}
    if mono:
        # result['mono'][lib][chrom][A] = count
        result['mono'] = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    if dinuc:
        result['dinuc'] = {}
        for i in dist:
            result['dinuc'][i] = {}
            result['dinuc'][i]['nr'] = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
            result['dinuc'][i]['rn'] = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    if trinuc:
        result['trinuc'] = {}
        result['trinuc']['nnr'] = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
        result['trinuc']['nrn'] = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
        result['trinuc']['rnn'] = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    # calculate cache length:
    cache_len = 0
    if trinuc:
        cache_len = 2
    if dinuc:
        cache_len = max(dist+[cache_len])
    # read
    cr = None
    genome = ''
    for l in gr:
        l = l.rstrip('\n')
        # header
        if l[0] == '>':
            # calculate chromosome
            calc_chrom(cr, genome, ribos, libs, mono, dinuc, trinuc, dist, result)
            # initialize
            cr = l[1:]
            genome = ''
        else:
            genome += l.rstrip('\n').upper()
    calc_chrom(cr, genome, ribos, libs, mono, dinuc, trinuc, dist, result)
    return result


# calculate for chromosome
def calc_chrom(cr, genome, ribos, libs, mono, dinuc, trinuc, dist, result):
    if not cr:
        return
    for pos, st, count in ribos[cr]:
        if mono:
            nuc = genome[pos-1]
            if st == '-':
                nuc = rc(nuc)
            add_ribo(nuc, libs, count, cr, result['mono'])
        if dinuc:
            for d in dist:
                if st == '+':
                    # nr
                    if pos - 1 - d >= 0:
                        nuc = genome[pos-1-d] + genome[pos-1]
                        add_ribo(nuc, libs, count, cr, result['dinuc'][d]['nr'])
                    # rn
                    if pos - 1 + d < len(genome):
                        nuc = genome[pos-1] + genome[pos-1 + d]
                        add_ribo(nuc, libs, count, cr, result['dinuc'][d]['rn'])
                else:
                    # nr
                    if pos - 1 + d < len(genome):
                        nuc = rc(genome[pos-1] + genome[pos-1+d])
                        add_ribo(nuc, libs, count, cr, result['dinuc'][d]['nr'])
                    # rn
                    if pos - 1 - d >= 0:
                        nuc = rc(genome[pos-1-d] + genome[pos-1])
                        add_ribo(nuc, libs, count, cr, result['dinuc'][d]['rn'])
        if trinuc:
            if st == '+':
                # nnr
                if pos >=3:
                    nuc = genome[pos-3:pos]
                    add_ribo(nuc, libs, count, cr, result['trinuc']['nnr'])
                # nrn
                if 2<= pos <= len(genome) -1:
                    nuc = genome[pos-2:pos +1]
                    add_ribo(nuc, libs, count, cr, result['trinuc']['nrn'])
                # rnn
                if  pos <= len(genome) -2:
                    nuc = genome[pos-1:pos +2]
                    add_ribo(nuc, libs, count, cr, result['trinuc']['rnn'])
            else:
                # nnr
                if  pos <= len(genome) -2:
                    nuc = rc(genome[pos-1:pos +2])
                    add_ribo(nuc, libs, count, cr, result['trinuc']['nnr'])
                # nrn
                if 2<= pos <= len(genome) -1:
                    nuc = rc(genome[pos-2:pos +1])
                    add_ribo(nuc, libs, count, cr, result['trinuc']['nrn'])
                # rnn
                if pos >=3:
                    nuc = rc(genome[pos-3:pos])
                    add_ribo(nuc, libs, count, cr, result['trinuc']['rnn'])
    del ribos[cr]
    print(f'{cr} finished! Length = {len(genome)}')


# add ribos for one entry
def add_ribo(nuc, libs, counts, cr, result):
    for i in range(len(counts)):
        if counts[i]:
            result[libs[i]][cr][nuc] += counts[i]


# output
def output(results, outputbase):
    # build order
    order = {'mono':[], 'nr':[], 'rn':[], 'nnr':[], 'nrn':[], 'rnn':[]}
    base = ['A','C','G','T']
    for i in base:
        order['mono'].append(i)
        for j in base:
            order['nr'].append(j+i)
            order['rn'].append(i+j)
            for k in base:
                order['nnr'].append(k+j+i)
                order['nrn'].append(j+i+k)
                order['rnn'].append(i+j+k)
    # output
    for k,v in results.items():
        if k == 'mono':
            for lib, v1 in v.items():
                with open(generate_outputname(outputbase, lib) + '.mono', 'w') as fw:
                    fw.write('\t'.join(['Sample'] + order['mono']) + '\n')
                    for cr, v2 in v1.items():
                        fw.write(cr + '\t' + '\t'.join(list(map(str, [v2[x] for x in order['mono']]))) + '\n')
        if k == 'dinuc':
            for d, v00 in v.items():
                for o, v0 in v00.items():
                    if o == 'nr':
                        for lib, v1 in v0.items():
                            with open(generate_outputname(outputbase, lib) +  f'.{k}_d{d}_nr', 'w') as fw:
                                fw.write('\t'.join(['Sample'] + order['nr']) + '\n')
                                for cr, v2 in v1.items():
                                    fw.write(cr + '\t' + '\t'.join(list(map(str, [v2[x] for x in order['nr']]))) + '\n')
                    else:
                        for lib, v1 in v0.items():
                            with open(generate_outputname(outputbase, lib) +  f'.{k}_d{d}_rn', 'w') as fw:
                                fw.write('\t'.join(['Sample'] + order['rn']) + '\n')
                                for cr, v2 in v1.items():
                                    fw.write(cr + '\t' + '\t'.join(list(map(str, [v2[x] for x in order['rn']]))) + '\n')
        if k == 'trinuc':
            for o, v0 in v.items():
                if o == 'nnr':
                    for lib, v1 in v0.items():
                        with open(generate_outputname(outputbase, lib) + '.trinuc_nnr', 'w') as fw:
                            fw.write('\t'.join(['Sample'] + order['nnr']) + '\n')
                            for cr, v2 in v1.items():
                                fw.write(cr + '\t' + '\t'.join(list(map(str, [v2[x] for x in order['nnr']]))) + '\n')
                elif o=='nrn':
                    for lib, v1 in v0.items():
                        with open(generate_outputname(outputbase, lib) + '.trinuc_nrn', 'w') as fw:
                            fw.write('\t'.join(['Sample'] + order['nrn']) + '\n')
                            for cr, v2 in v1.items():
                                fw.write(cr + '\t' + '\t'.join(list(map(str, [v2[x] for x in order['nrn']]))) + '\n')
                else:
                    for lib, v1 in v0.items():
                        with open(generate_outputname(outputbase, lib) + '.trinuc_rnn', 'w') as fw:
                            fw.write('\t'.join(['Sample'] + order['rnn']) + '\n')
                            for cr, v2 in v1.items():
                                fw.write(cr + '\t' + '\t'.join(list(map(str, [v2[x] for x in order['rnn']]))) + '\n')


# generate output name
def generate_outputname(outputbase, lib):
    if not outputbase:
        return lib
    return outputbase + '_' + lib










