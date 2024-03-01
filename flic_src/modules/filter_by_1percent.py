import numpy as np
import os

from collections import defaultdict
from operator import itemgetter


def read_isoform_files(path_to_dir, thr_1st, thr_2nd):
    d_of_isoforms = defaultdict(list)
    d_good_iso = {}

    for file in os.listdir(path_to_dir):
        with open(f'{path_to_dir}{file}') as rep:
            for line in rep:
                *line_l, cov = line.strip('\n').split('\t')
                cov = int(cov)
                d_of_isoforms['\t'.join(line_l)].append(cov)

    for key, val in d_of_isoforms.items():
        if sum([x >= thr_1st for x in val]) > 1 and sum([x >= thr_2nd for x in val]) > 0:
            d_good_iso[key] = np.mean(val)
    return d_good_iso


def read_filt_iso_file(fname):
    d_of_iso_with_gene_ids = {}
    with open(fname) as inf:
        for line in inf:
            *isoform, iso_id = line.strip('\n').split('\t')
            isoform = '\t'.join(isoform)
            d_of_iso_with_gene_ids[isoform] = iso_id
    return d_of_iso_with_gene_ids


def calc_expr_by_genes(iso_expr, d_of_iso_with_gene_ids):
    expr_by_genes = defaultdict(int)
    for iso, expr in iso_expr.items():
        if iso in d_of_iso_with_gene_ids.keys():
            gene_id = '.'.join(d_of_iso_with_gene_ids[iso].split('.')[:-1])
            expr_by_genes[gene_id] += expr
    return expr_by_genes


def filtering_iso(iso_expr, d_of_iso_with_gene_ids, expr_by_genes):
    s_good_iso = set()
    for iso, cur_iso_expr in iso_expr.items():
        if iso in d_of_iso_with_gene_ids.keys():
            gene_id = '.'.join(d_of_iso_with_gene_ids[iso].split('.')[:-1])
            cur_gene_expr = expr_by_genes[gene_id]
            if cur_iso_expr >= cur_gene_expr * 0.01:
                s_good_iso.add((iso, d_of_iso_with_gene_ids[iso]))
    return s_good_iso


def write_res(ouf_name, s_good_iso):
    with open(ouf_name, 'w') as ouf:
        s_good_iso = sorted(s_good_iso, key=itemgetter(1))
        for elem in s_good_iso:
            ouf.write(f'{elem[0]}\t{elem[1]}\n')


def filter_by_1percent(iso_dir, iso_thr1, iso_thr2, final_iso_file, out_dir):
    ouf_name = f'{out_dir}filt_1perc_{os.path.basename(final_iso_file)}'
    iso_expr = read_isoform_files(iso_dir, iso_thr2, iso_thr1)
    d_of_iso_with_gene_ids = read_filt_iso_file(final_iso_file)
    expr_by_genes = calc_expr_by_genes(iso_expr, d_of_iso_with_gene_ids)
    s_good_iso = filtering_iso(iso_expr, d_of_iso_with_gene_ids, expr_by_genes)
    write_res(ouf_name, s_good_iso)
