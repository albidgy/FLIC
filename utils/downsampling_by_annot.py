import concurrent.futures

import numpy as np
import os
import re
import shutil

from collections import defaultdict
from joblib import Parallel, delayed

D_OF_GOOD_FLAGS = {'0': '+', '16': '-'}


def read_annot_file(gtf_file):
    d_of_gene_coords = defaultdict(lambda: defaultdict(list))
    d_of_genes = defaultdict(lambda: defaultdict(str))

    with open(gtf_file) as annot:
        for line in annot:
            if line[0] == '#':
                continue

            line_l = line.strip('\n').split('\t')
            feat_type = line_l[2]

            if feat_type == 'gene':
                chrom_and_orientation = f'{line_l[0]}*{line_l[6]}'
                start_pos = int(line_l[3])
                end_pos = int(line_l[4])
                gene_id = re.findall(r'gene_id "(.+?)";', line_l[8])[0]

                for idx in range(start_pos, end_pos):
                    d_of_gene_coords[chrom_and_orientation][idx].append((start_pos, end_pos))
                d_of_genes[chrom_and_orientation][(start_pos, end_pos)] = gene_id
    return d_of_gene_coords, d_of_genes


def drop_multimapping_reads(filename, downsampling_outdir):
    d_of_reads = defaultdict(list)
    system_lines = []
    out_filename = f'{downsampling_outdir}unique_map_{os.path.basename(filename)}'
    with open(filename) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            if line[0] == '@':
                system_lines.append(line)
                continue
            read_id = line_l[0]
            assembly_type = line_l[1]
            d_of_reads[read_id].append((assembly_type, line))

    with open(out_filename, 'w') as ouf:
        ouf.write(''.join(system_lines))
        for key, val in d_of_reads.items():
            if len(val) > 1:
                continue
            assembly_type, line = val[0]
            if assembly_type not in D_OF_GOOD_FLAGS.keys():
                continue
            ouf.write(''.join(line))

    return out_filename


def split_cigar(cigar_str):
    l_of_splice_sites = []
    tmp_val = []
    for elem in cigar_str:
        if elem.isdigit():
            tmp_val.append(elem)
        else:
            l_of_splice_sites.append((elem, ''.join(tmp_val)))
            tmp_val = []
    return l_of_splice_sites


def get_correct_intron_and_end_pos(cigar_str, start_coord):
    counter_nucls = 0
    pos_of_introns = []
    l_of_splice_sites = split_cigar(cigar_str)

    for idx in range(len(l_of_splice_sites)):
        cigar, n_nucl = l_of_splice_sites[idx]
        n_nucl = int(n_nucl)
        if cigar == 'S':
            if idx == 0 or idx == len(l_of_splice_sites) - 1:
                continue
        elif cigar in 'I':
            continue
        elif cigar == 'N':
            pos_of_introns.append((start_coord + counter_nucls, start_coord + counter_nucls + n_nucl - 1))
        counter_nucls += n_nucl
    end_coord = start_coord + counter_nucls
    return end_coord, pos_of_introns


def check_intersection_with_gene(read_coords, genes_coords):
    d_of_cur_prop = {}
    start_read, end_read = read_coords
    for cur_broads in genes_coords:
        start_ref, end_ref = cur_broads
        max_start = max(start_ref, start_read)
        min_end = min(end_ref, end_read)
        min_length = min(end_ref - start_ref, end_read - start_read)
        prop_intersection = max(min_end - max_start, 0) / min_length

        if prop_intersection >= 0.8:
            d_of_cur_prop[prop_intersection] = cur_broads

    if d_of_cur_prop != {}:
        return d_of_cur_prop[max(d_of_cur_prop.keys())]
    return None


def make_sampling(pid, d_of_gene_reads_cov, tmp_outdir, max_n_reads, min_n_reads):
    s_downsampled_reads = []
    for key, val in d_of_gene_reads_cov.items():
        if len(val) < min_n_reads:
            continue
        elif len(val) <= max_n_reads:
            s_downsampled_reads.extend(val)
        else:
            s_downsampled_reads.extend(np.random.choice(val, max_n_reads, replace=False))

    with open(f'{tmp_outdir}{str(pid)}.txt', 'w') as ouf:
        ouf.write('\n'.join(s_downsampled_reads) + '\n')


def get_sampled_read_names(pid, sam_file, d_of_gene_coords, d_of_genes, tmp_outdir, max_n_reads, min_n_reads):
    d_of_gene_reads_cov = defaultdict(list)
    with open(sam_file) as sam:
        for line in sam:
            if line[0] == '@':
                continue
            line_l = line.strip('\n').split('\t')
            read_id = line_l[0]
            chrom = line_l[2]
            strand = D_OF_GOOD_FLAGS[line_l[1]]
            chrom_and_strand = f'{chrom}*{strand}'
            if chrom_and_strand not in d_of_gene_coords.keys():
                continue
            start_read = int(line_l[3])
            cigar_str = line_l[5]
            end_read, _ = get_correct_intron_and_end_pos(cigar_str, int(line_l[3]))
            middle_coord_of_read = int((end_read - start_read) / 2 + start_read)

            if middle_coord_of_read in d_of_gene_coords[chrom_and_strand].keys():
                gene_broads = check_intersection_with_gene((start_read, end_read),
                                                           d_of_gene_coords[chrom_and_strand][middle_coord_of_read])
                if gene_broads is not None:
                    d_of_gene_reads_cov[d_of_genes[chrom_and_strand][gene_broads]].append(read_id)

    make_sampling(pid, d_of_gene_reads_cov, tmp_outdir, max_n_reads, min_n_reads)


def prepare_data_for_multiprocessing(d_of_gene_coords, d_of_genes, n_proc):
    d_of_gene_coords_by_proc = defaultdict(dict)
    d_of_genes_by_proc = defaultdict(dict)
    counter = 0

    if n_proc > len(d_of_gene_coords):
        n_proc = len(d_of_gene_coords)

    for key, val in d_of_gene_coords.items():
        d_of_gene_coords_by_proc[counter % n_proc][key] = d_of_gene_coords[key]
        d_of_genes_by_proc[counter % n_proc][key] = d_of_genes[key]
        counter += 1
    return d_of_gene_coords_by_proc, d_of_genes_by_proc


def filter_sam_file(uniq_sam, tmp_outdir):
    s_downsampled_reads = set()
    final_outname = uniq_sam.replace('unique_map_', '')

    for file in os.listdir(tmp_outdir):
        with open(f'{tmp_outdir}{file}') as inf:
            for line in inf:
                s_downsampled_reads.add(line.strip('\n'))

    with open(uniq_sam) as sam_file:
        for line in sam_file:
            if line[0] == '@':
                with open(final_outname, 'a') as ouf:
                    ouf.write(line)
            else:
                read_id = line.split('\t')[0]
                if read_id in s_downsampled_reads:
                    with open(final_outname, 'a') as ouf:
                        ouf.write(line)


def run_downsampling(need_make_downsampling, sam_file, gtf_file, num_threads, max_n_reads, min_n_reads, common_outdir):
    downsampling_outdir = common_outdir + 'downsampled/'
    if not os.path.exists(downsampling_outdir):
        os.mkdir(downsampling_outdir)

    d_of_gene_coords, d_of_genes = read_annot_file(gtf_file)
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:  # need for clearing memory
        uniq_map_sam = executor.submit(drop_multimapping_reads, sam_file, downsampling_outdir,).result()

    if need_make_downsampling:
        tmp_outdir = downsampling_outdir + 'tmp/'
        os.mkdir(tmp_outdir)
        d_of_gene_coords_by_proc, d_of_genes_by_proc = prepare_data_for_multiprocessing(d_of_gene_coords,
                                                                                        d_of_genes, num_threads)
        _ = Parallel(n_jobs=num_threads)(delayed(get_sampled_read_names)(pid, uniq_map_sam,
                                                                         d_of_gene_coords_by_proc[pid],
                                                                         d_of_genes_by_proc[pid], tmp_outdir,
                                                                         max_n_reads, min_n_reads)
                                         for pid in d_of_gene_coords_by_proc.keys())

        filter_sam_file(f'{uniq_map_sam}', tmp_outdir)
        shutil.rmtree(tmp_outdir)
    return downsampling_outdir, uniq_map_sam
