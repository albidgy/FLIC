import networkx as nx
import os
import re
import shutil

from collections import defaultdict
from joblib import Parallel, delayed


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
            d_good_iso[key] = len(val)
    return d_good_iso


def make_d_borders_by_iso(d_of_isoforms):
    d_border_iso_with_splice_sites = {}

    for key in d_of_isoforms.keys():
        line_l = key.split('\t')
        chrom_and_orientation = f'{line_l[0]}*{line_l[1]}'
        start = int(line_l[2].split('-')[0])
        stop = int(line_l[4].split('-')[1])
        splice_sites = line_l[3].split(';')
        if splice_sites == ['']:
            splice_sites = []

        if chrom_and_orientation not in d_border_iso_with_splice_sites.keys():
            d_border_iso_with_splice_sites[chrom_and_orientation] = {}
        d_border_iso_with_splice_sites[chrom_and_orientation][(start, stop)] = splice_sites
    return d_border_iso_with_splice_sites


def calc_similaritivity(l_of_borders):
    l_of_pair_similaritivity = []

    for idx in range(len(l_of_borders)):
        for jdx in range(len(l_of_borders)):
            max_start = max(l_of_borders[idx][0], l_of_borders[jdx][0])
            min_end = min(l_of_borders[idx][1], l_of_borders[jdx][1])
            min_length = min(l_of_borders[idx][1] - l_of_borders[idx][0],
                             l_of_borders[jdx][1] - l_of_borders[jdx][0])

            prop_intersection = max(min_end - max_start, 0) / min_length
            l_of_pair_similaritivity.append((l_of_borders[idx], l_of_borders[jdx], round(prop_intersection, 2)))
    return l_of_pair_similaritivity


def create_graph(l_of_borders, l_of_pair_similaritivity):
    graph = nx.Graph()
    graph.add_nodes_from(l_of_borders)
    graph.add_weighted_edges_from(l_of_pair_similaritivity)
    bad_scores = list(filter(lambda e: e[2] < 0.8, (e for e in graph.edges.data('weight'))))
    bad_ids = list(e[:2] for e in bad_scores)
    graph.remove_edges_from(bad_ids)
    return graph


def extract_longest_completted_subgraphs(graph):
    result = []
    len_longest = 100
    while len_longest > 1:
        longest_sub = []
        len_longest = 100
        l_of_complete_subgraphs = list(nx.find_cliques(graph))
        if not l_of_complete_subgraphs:
            len_longest = 0
        for sub in l_of_complete_subgraphs:
            if len_longest == 100:
                len_longest = len(sub)
                longest_sub = sub
            if len(sub) > len_longest:
                longest_sub = sub
                len_longest = len(sub)
        if len_longest > 1:
            result.append(longest_sub)
        graph.remove_nodes_from(longest_sub)
    result.extend(l_of_complete_subgraphs)
    return result


def find_gene_borders_and_splice_sites(result_itteration, d_cur_borders):
    d_new_borders = {}
    for merged_coords in result_itteration:
        new_start_coord = min(merged_coords, key=lambda x: x[0])[0]
        new_stop_coord = max(merged_coords, key=lambda x: x[1])[1]
        d_new_borders[(new_start_coord, new_stop_coord)] = []

        for coord in merged_coords:
            cur_start, cur_stop = coord
            cur_splice_sites = d_cur_borders[(cur_start, cur_stop)]
            d_new_borders[(new_start_coord, new_stop_coord)].extend(cur_splice_sites)

    for key, val in d_new_borders.items():
        d_new_borders[key] = list(set(val))
    return d_new_borders


def run_clustering_by_graph(d_cur_borders):
    l_borders = list(d_cur_borders.keys())

    l_of_pair_similaritivity = calc_similaritivity(l_borders)
    graph = create_graph(l_borders, l_of_pair_similaritivity)
    result_itteration = extract_longest_completted_subgraphs(graph)
    d_new_borders = find_gene_borders_and_splice_sites(result_itteration, d_cur_borders)
    return d_new_borders


def write_genes(d_of_gene_borders, chrom_and_orientation, output_name):
    chrom, orientation = chrom_and_orientation.split('*')
    with open(output_name, 'a') as ouf:
        for key, val in d_of_gene_borders.items():
            gene_start, gene_stop = key
            splice_sites_str = ';'.join(val)
            ouf.write(f'{chrom}\t{orientation}\t{str(gene_start)}\t{splice_sites_str}\t{str(gene_stop)}\n')


def make_genes(pid, d_border_iso_with_splice_sites, ouf_filename):
    out_filename = f'{ouf_filename}_{str(pid)}.tsv'
    for key, val in d_border_iso_with_splice_sites.items():
        prev_l_of_borders = []
        d_cur_borders = d_border_iso_with_splice_sites[key]

        while prev_l_of_borders != list(d_cur_borders.keys()):
            prev_l_of_borders = list(d_cur_borders.keys())
            d_cur_borders = dict(sorted(run_clustering_by_graph(d_cur_borders).items()))
        write_genes(d_cur_borders, key, out_filename)


def prepare_data_for_multiprocessing(d_border_iso_with_splice_sites, n_proc):
    data_by_proc = {}
    counter = 0

    if n_proc > len(d_border_iso_with_splice_sites):
        n_proc = len(d_border_iso_with_splice_sites)
    for proc in range(n_proc):
        data_by_proc[proc] = {}

    for key, val in d_border_iso_with_splice_sites.items():
        data_by_proc[counter % n_proc][key] = d_border_iso_with_splice_sites[key]
        counter += 1
    return data_by_proc


def concate_genes(out_filename, path_to_tmp_dir):
    for file in os.listdir(path_to_tmp_dir):
        with open(f'{path_to_tmp_dir}{file}') as inf:
            for line in inf:
                with open(out_filename, 'a') as ouf:
                    ouf.write(line)


def create_genes(iso_dir, num_threads, thr_1st, thr_2nd, output_dir):
    genes_dir = output_dir + 'genes/'
    os.mkdir(genes_dir)
    tmp_outdir = genes_dir + 'tmp/'
    os.mkdir(tmp_outdir)
    out_filename = 'genes'

    d_good_iso = read_isoform_files(iso_dir, thr_1st, thr_2nd)
    d_border_iso_with_splice_sites = make_d_borders_by_iso(d_good_iso)
    data_by_proc = prepare_data_for_multiprocessing(d_border_iso_with_splice_sites, num_threads)

    _ = Parallel(n_jobs=num_threads)(delayed(make_genes)(pid, d_border_iso_splitted, tmp_outdir + out_filename)
                                     for pid, d_border_iso_splitted in data_by_proc.items())
    concate_genes(f'{genes_dir}{out_filename}.tsv', tmp_outdir)
    shutil.rmtree(tmp_outdir)
    return genes_dir


def read_annot_file(annot_file):
    d_of_annot_genes = defaultdict(dict)
    with open(annot_file) as inf:
        for line in inf:
            if line[0] == '#':
                continue
            line_l = line.strip('\n').split('\t')
            feat_id = line_l[2]
            if feat_id == 'gene':
                chrom_and_orientation = f'{line_l[0]}*{line_l[6]}'
                start = int(line_l[3])
                stop = int(line_l[4])
                gene_id = re.findall(r'gene_id "(.+?)";', line_l[8])[0]

                d_of_annot_genes[chrom_and_orientation][(start, stop)] = gene_id
    return d_of_annot_genes


def intersection_with_genes(gene_coords, d_of_ref):
    start_gene, stop_gene = gene_coords

    for key, gene_id in d_of_ref.items():
        start_ref, stop_ref = key

        max_start = max(start_ref, start_gene)
        min_end = min(stop_ref, stop_gene)
        min_length = min(stop_ref - start_ref,
                         stop_gene - start_gene)
        prop_intersection = max(min_end - max_start, 0) / min_length

        if prop_intersection >= 0.8:
            return gene_id
    return 'unassigned_gene'


def move_geneids_from_annot(annot_file, genes_file, out_dir):
    out_fname = f'{out_dir}final_results/{os.path.basename(genes_file)}'
    os.mkdir(f'{out_dir}final_results/')

    d_annot_file = read_annot_file(annot_file)
    with open(out_fname, 'w') as ouf:
        ouf.write('')

    with open(genes_file) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            chrom_and_orientation = f'{line_l[0]}*{line_l[1]}'
            start = int(line_l[2])
            stop = int(line_l[4])

            gene_id = intersection_with_genes((start, stop), d_annot_file[chrom_and_orientation])
            line_l.append(gene_id)

            with open(out_fname, 'a') as ouf:
                ouf.write('\t'.join(line_l) + '\n')

    return out_fname
