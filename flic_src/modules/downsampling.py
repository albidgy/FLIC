import concurrent.futures
import logging
import os
import re
import shutil
import sys
from collections import defaultdict

import numpy as np
from joblib import Parallel, delayed

D_OF_GOOD_FLAGS = {'0': '+', '16': '-'}


def read_annot_file(gtf_file):
    d_of_gene_coords = defaultdict(lambda: defaultdict(list))
    d_of_genes = defaultdict(lambda: defaultdict(str))
    is_find_feat_gene = False

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
                gene_id = re.findall(r'gene_id "(.+?)"', line_l[8])[0]
                is_find_feat_gene = True

                for idx in range(start_pos, end_pos):
                    d_of_gene_coords[chrom_and_orientation][idx].append((start_pos, end_pos))
                d_of_genes[chrom_and_orientation][(start_pos, end_pos)] = gene_id

    if not is_find_feat_gene:
        print('Error: The annotation needs to contain information about the gene structure (3rd column: feature). '
              'See the GTF file in the example directory',
              file=sys.stderr)
        sys.exit(1)
    return d_of_gene_coords, d_of_genes


def make_d_of_chroms_and_d_of_read_coords_wo_annot(sam_file):
    d_of_chroms = defaultdict(lambda: np.zeros(chrom_len, dtype=int))
    d_of_reads_coords = defaultdict(lambda: defaultdict(dict))
    d_of_chrom_lens = {}

    with open(sam_file) as map_sam:
        for line in map_sam:
            line_l = line.strip('\n').split('\t')
            if line[0] == '@':
                if line.startswith('@SQ'):
                    chromosome_name = line_l[1].split('SN:')[1]
                    chromosome_len = int(line_l[2].split('LN:')[1]) + 3  # needed for correct work with extranuclear DNA
                    d_of_chrom_lens[chromosome_name] = chromosome_len
                else:
                    continue
            else:
                read_id = line_l[0]
                chrom = line_l[2]
                strand = D_OF_GOOD_FLAGS[line_l[1]]
                chrom_and_strand = f'{chrom}_{strand}'
                chrom_len = d_of_chrom_lens[chrom]
                start_ref = int(line_l[3]) + 1  # needed for correct work with extranuclear DNA
                cigar_str = line_l[5]
                end_ref, _ = get_correct_intron_and_end_pos(cigar_str, int(line_l[3]))

                d_of_chroms[chrom_and_strand][start_ref:end_ref + 2] += 1  # needed for correct work with extranuclear DNA

                split_starts_by_100_000 = start_ref // 100_000
                if split_starts_by_100_000 not in d_of_reads_coords[chrom_and_strand].keys():
                    d_of_reads_coords[chrom_and_strand][split_starts_by_100_000] = {}
                d_of_reads_coords[chrom_and_strand][split_starts_by_100_000][read_id] = (start_ref, end_ref)
    return d_of_chroms, d_of_reads_coords


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
        elif cigar == 'I':
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


def add_reads_into_black_list(idx, d_of_reads_coords, black_list, max_n_reads):
    l_of_searched_reads = []
    ref_coord = idx // 100_000
    l_of_keys = [ref_coord - 1, ref_coord, ref_coord + 1]

    for cur_itter in l_of_keys:
        if cur_itter in d_of_reads_coords.keys():
            for read_id, val in d_of_reads_coords[cur_itter].items():
                if read_id not in black_list and val[0] <= idx <= val[1]:
                    l_of_searched_reads.append(read_id)
    if len(l_of_searched_reads) > max_n_reads:
        black_list.update(np.random.choice(l_of_searched_reads, len(l_of_searched_reads) - max_n_reads, replace=False))
    return black_list


def cut_cov_by_max_thr_wo_annot(pid, d_of_chroms, d_of_reads_coords, tmp_outdir, max_n_reads):
    black_list = set()
    for key, l_of_genome_cov in d_of_chroms.items():
        for idx in range(len(l_of_genome_cov)):
            if l_of_genome_cov[idx] > max_n_reads:
                black_list = add_reads_into_black_list(idx, d_of_reads_coords[key], black_list, max_n_reads)

    with open(f'{tmp_outdir}{pid}.txt', 'w') as ouf:
        ouf.write('\n'.join(black_list) + '\n')


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


def filter_sam_file_wo_annot(uniq_sam, tmp_outdir):
    s_black_list = set()
    final_outname = uniq_sam.replace('unique_map_', '')

    for file in os.listdir(tmp_outdir):
        with open(f'{tmp_outdir}{file}') as inf:
            for line in inf:
                s_black_list.add(line.strip('\n'))

    with open(uniq_sam) as sam_file:
        for line in sam_file:
            if line[0] == '@':
                with open(final_outname, 'a') as ouf:
                    ouf.write(line)
            else:
                read_id = line.split('\t')[0]
                if read_id not in s_black_list:
                    with open(final_outname, 'a') as ouf:
                        ouf.write(line)


def run_downsampling(need_make_downsampling, ref_annot, sam_file, gtf_file,
                     num_threads, max_n_reads, min_n_reads, common_outdir):
    logging.info(f'    Downsample reads by max threshold {str(max_n_reads)} reads per gene')
    downsampling_outdir = common_outdir + 'downsampled/'
    if not os.path.exists(downsampling_outdir):
        os.mkdir(downsampling_outdir)

    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:  # need for clearing memory
        uniq_map_sam = executor.submit(drop_multimapping_reads, sam_file, downsampling_outdir, ).result()

    if need_make_downsampling:
        tmp_outdir = downsampling_outdir + 'tmp/'
        os.mkdir(tmp_outdir)
        if ref_annot is not None:
            d_of_gene_coords, d_of_genes = read_annot_file(gtf_file)
            d_of_gene_coords_by_proc, d_of_genes_by_proc = prepare_data_for_multiprocessing(d_of_gene_coords,
                                                                                            d_of_genes, num_threads)
            _ = Parallel(n_jobs=num_threads)(delayed(get_sampled_read_names)(pid, uniq_map_sam,
                                                                             d_of_gene_coords_by_proc[pid],
                                                                             d_of_genes_by_proc[pid], tmp_outdir,
                                                                             max_n_reads, min_n_reads)
                                             for pid in d_of_gene_coords_by_proc.keys())

            filter_sam_file(f'{uniq_map_sam}', tmp_outdir)
        else:
            d_of_chroms, d_of_reads_coords = make_d_of_chroms_and_d_of_read_coords_wo_annot(uniq_map_sam)
            d_of_chroms_by_proc, d_of_reads_coords_by_proc = prepare_data_for_multiprocessing(d_of_chroms,
                                                                                              d_of_reads_coords,
                                                                                              num_threads)
            _ = Parallel(n_jobs=num_threads)(delayed(cut_cov_by_max_thr_wo_annot)(pid,
                                                                                  d_of_chroms_by_proc[pid],
                                                                                  d_of_reads_coords_by_proc[pid],
                                                                                  tmp_outdir,
                                                                                  max_n_reads)
                                             for pid in d_of_chroms_by_proc.keys())

            filter_sam_file_wo_annot(uniq_map_sam, tmp_outdir)

        shutil.rmtree(tmp_outdir)
    return downsampling_outdir, uniq_map_sam
