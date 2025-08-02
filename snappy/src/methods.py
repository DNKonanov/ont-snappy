
from Bio.SeqIO import parse
from snappy.src.seq_processing import gen_variants, letter_codes, letter_codes_rev
import re    
from scipy.stats import mode
import pandas as pd
from itertools import product
import polars as pl
from tqdm import tqdm
import numpy as np
from scipy.stats import chi2_contingency

from Bio.Seq import reverse_complement

import warnings
warnings.filterwarnings('ignore')

MODKIT_COLNAMES = (
        'contig',
        'start',
        'end',
        'mod_type',
        'score',
        'strand',
        'start_pos',
        'end_pos',
        'color',
        'n_valid',
        'frac_mod',
        'n_mod',
        'n_canon',
        'n_other',
        'n_delete',
        'n_fail',
        'n_diff',
        'n_nocall',
    )


canonical_bases = ['A', 'C', 'G', 'T']

re_to_degenerate = {
    '[A,]': 'A',
    '[C,]': 'C',
    '[G,]': 'G',
    '[T,]': 'T',
    '[A,C,]': 'M',
    '[A,G,]': 'R',
    '[A,T,]': 'W',
    '[C,G,]': 'S',
    '[C,T,]': 'Y',
    '[G,T,]': 'K',
    '[A,C,G,]': 'V',
    '[A,C,T,]': 'H',
    '[A,G,T,]': 'D',
    '[C,G,T,]': 'B',
    '[A,C,G,T,]': 'N',
}

degenerate_weght = {
    'A': 1,
    'C': 1,
    'G': 1,
    'T': 1,
    'M': 2,
    'R': 2,
    'W': 2,
    'S': 2,
    'Y': 2,
    'K': 2,
    'V': 3,
    'H': 3,
    'D': 3,
    'B': 3,
    'N': 4,
    '.': 4,
}



degenerate_to_re  = {
    re_to_degenerate[k]: k for k in re_to_degenerate
}

degenerate_bases = list(degenerate_to_re.keys())

OFFSETS_RANGE = list(range(0,12))

def add_contexts(df, reference):
    print('Adding contexts...')
    contexts = []
    for i in tqdm(range(len(df))):
        
        contig, start = str(df['contig'][i]), df['start'][i]
        if df['strand'][i] == '+':
            contexts.append(reference[contig][start - 15:start + 16])
        elif df['strand'][i]== '-':
            contexts.append(
                reverse_complement(reference[contig][start - 15:start + 16])
            )


    df = df.with_columns(
        pl.Series(name="context", values=contexts),
    )
    return df



def unite_motifs(m_list):    
    m_list = np.array([[s for s in m] for m in m_list]).T
    
    keys = [tuple(sorted(list(set(row)))) for row in m_list]
    
    new_motif = ''
    for k in keys:
        try:
            new_motif += letter_codes_rev[k]
        except KeyError:
            new_motif += k[0]
    return new_motif


def gen_kmers(repeat=3):
    
    kmers = [''.join(m) for m in product(canonical_bases, repeat=repeat)]
    return kmers 

kmers = gen_kmers()

def get_kmer_pair_keys(kmers):
    
    keys = list(product(kmers, repeat=2))
    return keys

def parse_modkit(file,):

    print('Read modkit table...')

    modkit_df = pl.read_csv(
        file, 
        separator = '\t', 
        has_header = False,
        new_columns = MODKIT_COLNAMES,
    )

    return modkit_df


def parse_fasta(file):
    fasta = parse(file, 'fasta')
    reference = {}
    for rec in fasta:
        reference[rec.description.split(' ')[0]] = str(rec.seq)
    return reference
    


def construct_control_graph(reference):

    kmers = gen_kmers()


    dict_freq_cntr = {}
    edge_count_cntr = 0

    print('Constructing control graph...')
    for offset in tqdm(OFFSETS_RANGE):
        
        dict_freq_cntr[offset] = {
            kmer_key: 0 for kmer_key in get_kmer_pair_keys(kmers)
        }
        
        
        for contig in reference.keys():
            for ref in (reference[contig], reverse_complement(reference[contig])):
                for i in range(len(ref)-4-offset):

                    source = ref[i:i+3]
                    target = ref[i+1+offset:i+4+offset]
                    try:
                        dict_freq_cntr[offset][(source, target)] += 1
                        edge_count_cntr += 1
                    except KeyError:
                        continue

    
    return dict_freq_cntr, edge_count_cntr


def get_dict_freq_exp(df, OFFSETS_RANGE):


    dict_freq_exp = {}
    edge_count_exp = 0

    print('Constructing exp graph...')
    for offset in tqdm(OFFSETS_RANGE):

        dict_freq_exp[offset] = {
            kmer_key: 0 for kmer_key in get_kmer_pair_keys(kmers)
        }

        for read in df['context']:

            for i in range(0 + offset, len(read) - 4 - offset):

                source = read[i:i+3]
                target = read[i+1+offset:i+4+offset]
                try:
                    dict_freq_exp[offset][(source, target)] += 1
                    edge_count_exp += 1
                except KeyError:
                    continue
                    
    return dict_freq_exp


def get_chi_stats(dict_freq_exp, dict_freq_cntr, min_freq_ratio = 2):
    

    chi_stat_tables = {}
    kmers = gen_kmers()
    
    print(f'Search for differential edges...')
    for _offset in tqdm(OFFSETS_RANGE):
        chi_stat_table = pd.DataFrame(index=kmers, columns=kmers)

        edge_count_exp = np.sum(list(dict_freq_exp[_offset].values()))
        edge_count_cntr = np.sum(list(dict_freq_cntr[_offset].values()))

        for k1 in chi_stat_table.columns:
            for k2 in chi_stat_table.index:
                
                f1, s1, f2, s2 = dict_freq_exp[_offset][(k1, k2)], edge_count_exp, dict_freq_cntr[_offset][(k1, k2)], edge_count_cntr

                

                if (f1/s1) / (f2/s2) < min_freq_ratio:
                    statistics = 0
                
                else:
                    _2x2_table = [
                        [f1, s1],
                        [f2, s2],
                    ]

                    try:
                        statistics = chi2_contingency(
                            _2x2_table
                        )[0]
                    except ValueError:
                        statistics = 0

                chi_stat_table[k1][k2] = statistics

        chi_stat_tables[_offset] = chi_stat_table
    return chi_stat_tables
        


def get_initial_max_edge(chi_stat_tables):


    max_edge = None
    max_confidence = -1

    for offset in OFFSETS_RANGE[:]:
        chi_stat_table = chi_stat_tables[offset]

        for k1 in chi_stat_table.columns:
            for k2 in chi_stat_table.index:
                confidence = chi_stat_table[k1][k2]
                if confidence > max_confidence:

                    max_edge = ((k1, k2), confidence, offset)
                    max_confidence = confidence

    return max_edge




def get_re_template(k1, k2,offset):
    
    if offset == 0:
        re_template = k1 + k2[-1:]
    elif offset == 1:
        re_template = k1 + k2[-2:]
    else:
        re_template = k1 + '.' * int((offset - len(k1)) + 1) + k2
        
    return re_template


def get_best_positions(pattern, contexts):
    
    p = re.compile(pattern)
    starts = []
    ends = []
    
    for context in contexts:
        
        for m in p.finditer(context):
            starts.append(m.start())
            ends.append(m.end())
    
    return mode(starts)[0], mode(ends)[0]



def get_modality(sdf, template, upper_thrs = 40, lower_thrs = 10, metric_thrs=0.87):
    
    t_variants = gen_variants(template)
    dists = []
    lens = []
    
    for v in t_variants:
        ssdf = sdf.filter(sdf['context'].str.contains(v))
        dists.append(list(ssdf['frac_mod']))
        lens.append(len(ssdf))
    
    max_repr = max(lens)
    dist = []
    for i, d in enumerate(dists):
        if lens[i] == 0:
            return False, 0.0
        dist += d * int(np.round(max_repr/lens[i]))

    dist = np.array(dist)
    
    if len(dist) == 0:
        return False, 0.0
    
    upper = len(dist[dist >= upper_thrs])
    lower = len(dist[dist <= lower_thrs])
    
    try:
        metric = upper/(upper + lower)
    except ZeroDivisionError:
        return False, 0.0
    
    if metric > metric_thrs:
        return True, metric
    
    else:
        return False, metric
    

def get_adjacent_positions(template):
    
    adjacent_positions = set([])
    for i, base in enumerate(template):
        if i == 0 or i == len(template) - 1:
            continue
        
        if base != '.': continue
        if template[i+1] != '.' or template[i-1] != '.':
            adjacent_positions.add(i)
            
    return sorted(list(adjacent_positions))
    
def get_non_N_positions(template):
    non_N_positions = []
    
    for i, base in enumerate(template):
        if base != '.':
            non_N_positions.append(i)
            
    return non_N_positions


def get_re_from_degenerate(template):
    
    new_template = ''
    for i, base in enumerate(template):
        if base == '.' or base in canonical_bases:
            new_template += base
        else:
            new_template += degenerate_to_re[base]
    return new_template        
    
def get_degenerate_from_re(re_template):
    
    new_template = re_template
    
    for key in re_to_degenerate:
        new_template = new_template.replace(key, re_to_degenerate[key])
    return new_template
    
    





def gen_position_specific_variants(template, pos):

    symbol = template[pos]

    s_templates = []
    for t in letter_codes[symbol]:
        s_templates.append(
            template[:pos] + t + template[pos+1:]
        )
    return s_templates


def is_degenerate(base):
    if base == '.' or base == 'N': return True
    return False

def curate_1_links(template, df,):

    positions_to_curate = []
    for i, base in enumerate(template):

        if  i == 0 or i == 31:
            continue

        if is_degenerate(base) and template[i-1] != '.' and template[i+1] !='.':
            positions_to_curate.append(i)

    if len(positions_to_curate) == 0:
        return template

    sdf = df.filter(df['context'].str.contains(template))
    
    for p in positions_to_curate:
        bases = tuple(sorted(list(set([seq[p] for seq in sdf['context']]))))
        new_base = letter_codes_rev[bases]


        template = template[:p] + new_base + template[p+1:]

    return template





def curate_motif(template, central_base, modkit_df, metric_thrs=0.87):
    

    if template[15] != central_base:
        template = template[:15] + central_base + template[16:]

    _modkit_local = modkit_df.clone()
    sdf = _modkit_local.filter(_modkit_local['context'].str.contains(template))
    is_uniformic, metric = get_modality(sdf, get_degenerate_from_re(template))
    print(metric)
    
    print('Initial template:', template)
    
    # extend template if required
    while is_uniformic == False:
        
        adjacent_positions = get_adjacent_positions(template)
        non_N_positions = get_non_N_positions(template)
        
        temp_sdf = _modkit_local.filter(_modkit_local['context'].str.contains(template))
        
        hits = []
        
        for pos in adjacent_positions:
            for s in degenerate_bases:
                
                curated_template = template[:pos] + s + template[pos+1:]
                re_curated_template = get_re_from_degenerate(curated_template)
                
                sdf = temp_sdf.filter(temp_sdf['context'].str.contains(re_curated_template))
                is_uniformic, metric = get_modality(sdf, curated_template)

                hits.append(
                    (curated_template, sum([degenerate_weght[base] for base in curated_template]), metric, is_uniformic)
                )

        if True in [p[3] for p in hits]:
            
            hits = [p for p in hits if p[3] == True]
            template = sorted(hits, reverse=True, key=lambda x: x[2])[0][0]
            is_uniformic = True
            
        else:
            template = sorted(hits, reverse=True, key=lambda x: x[2])[0][0]
        
    print('Uniformic template:', template)
    print()
    
    
    # iteratively check non-N positions
    print('Iterative correction...')
    non_N_positions = get_non_N_positions(template)
    for pos in non_N_positions:

        if pos == 15:
            continue
        pos_hits = []

        for s in degenerate_bases:
            if s in canonical_bases:
                continue

            curated_template = template[:pos] + s + template[pos+1:]
            re_curated_template = get_re_from_degenerate(curated_template)

            sdf = _modkit_local.filter(_modkit_local['context'].str.contains(re_curated_template))
            is_uniformic, metric = get_modality(sdf, curated_template)

            print(curated_template, re_curated_template, is_uniformic, np.round(metric, 3))

            if is_uniformic:
                pos_hits.append(
                    (curated_template, sum([degenerate_weght[base] for base in curated_template]))
                )

        if len(pos_hits) > 0:
            template = sorted(pos_hits, reverse=True, key=lambda x: x[1])[0][0]
            subvariant_hits = []


            for t in gen_position_specific_variants(template, pos):
                t_sdf = _modkit_local.filter(_modkit_local['context'].str.contains(get_re_from_degenerate(t)))
                is_uniformic, metric = get_modality(t_sdf, t)

                if metric >= metric_thrs:
                    subvariant_hits.append(t)
            template = unite_motifs(subvariant_hits)


        corrected_segment = template[:pos + 1] + '.'*(31 - pos - 1)
        print(f'Local filtering for position {pos + 1}...')
        _modkit_local = _modkit_local.filter(_modkit_local['context'].str.contains(get_re_from_degenerate(corrected_segment)))

        

    
    print('Corrected motif:', template)
    
    return template
