import functools
import logging
import os

from joblib import Parallel, delayed


def get_d_illumina_splice_sites(illumina_file):
    d_of_illumina_sites = {}
    with open(illumina_file) as illumina_sites:
        for line in illumina_sites:
            line_l = line.strip('\n').split('\t')
            chrom = line_l[0]
            orientation = line_l[1]
            start = int(line_l[2])
            stop = int(line_l[3])

            if f'{chrom}*{orientation}' not in d_of_illumina_sites.keys():
                d_of_illumina_sites[f'{chrom}*{orientation}'] = {}
            split_by_starts_mln = start // 10_000
            if split_by_starts_mln not in d_of_illumina_sites[f'{chrom}*{orientation}'].keys():
                d_of_illumina_sites[f'{chrom}*{orientation}'][split_by_starts_mln] = []
            d_of_illumina_sites[f'{chrom}*{orientation}'][split_by_starts_mln].append((start, stop))
    return d_of_illumina_sites


def get_orientation(sam_flag):
    if sam_flag == 0:
        return '+'
    else:
        if bin(sam_flag)[-5] == '0':
            return '+'
        else:
            return '-'


def split_cigar(cigar_str):
    l_of_splice_sites = []
    tmp_val = []
    for elem in cigar_str:
        if elem in {'1', '2', '3', '4', '5', '6', '7', '8', '9', '0'}:
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
    end_coord = start_coord + counter_nucls - 1
    return end_coord, pos_of_introns


def find_potential_ill_sites(l_of_ill_sites, start, stop):
    l_of_site_intersecrions = []
    for ill_site in l_of_ill_sites:
        start_ill, stop_ill = ill_site
        start_ill_range = range(start_ill - 10, start_ill + 11)
        stop_ill_range = range(stop_ill - 10, stop_ill + 11)
        if start_ill + 11 <= stop:
            if start in start_ill_range and stop in stop_ill_range:
                l_of_site_intersecrions.append((start_ill, stop_ill))
        else:
            break
    return l_of_site_intersecrions


def find_closest_sites(pos_of_introns, d_of_illumina_sites_chr):
    res_l_correct_sites = []
    for cur_intron in pos_of_introns:
        closest_start_distance = 1000
        best_positions = tuple()
        start, stop = cur_intron
        l_of_potential_sites = []
        split_by_starts = start // 10_000

        for split_start in [split_by_starts - 1, split_by_starts, split_by_starts + 1]:
            if split_start in d_of_illumina_sites_chr.keys():
                l_of_potential_sites.extend(find_potential_ill_sites(sorted(d_of_illumina_sites_chr[split_start]),
                                                                     start, stop))

        if l_of_potential_sites:
            for intersected_site in l_of_potential_sites:
                start_intersected, stop_intersected = intersected_site
                cur_distance = abs(start - start_intersected) + abs(stop - stop_intersected)
                if cur_distance < closest_start_distance:
                    closest_start_distance = cur_distance
                    best_positions = intersected_site
            res_l_correct_sites.append(best_positions)

    return res_l_correct_sites


def read_sam_file(sam_file):
    sam_lines = []
    with open(sam_file) as filtered_sam:
        for line in filtered_sam:
            if line[0] == '@':
                continue
            line_l = line.strip('\n').split('\t')
            assembly_type = line_l[16][-1]
            if assembly_type not in ['P', 'I']:  # for checking
                continue

            read_id = line_l[0]
            orientation = get_orientation(int(line_l[1]))
            chrom = line_l[2]
            start_coord = int(line_l[3])
            cigar_str = line_l[5]
            sam_lines.append([read_id, orientation, chrom, start_coord, cigar_str])
    return sam_lines


def prepare_data_for_multiprocessing(sam_lines, n_proc):
    data_by_proc = {}
    counter = 0

    for proc in range(n_proc):
        data_by_proc[proc] = []

    for line_l in sam_lines:
        data_by_proc[counter % n_proc].append(line_l)
        counter += 1
    return data_by_proc


def extract_splice_sites(sam_lines, d_of_illumina_sites):
    l_res_lines = []
    for line_l in sam_lines:
        read_id, orientation, chrom, start_coord, cigar_str = line_l
        end_coord, pos_of_introns = get_correct_intron_and_end_pos(cigar_str, start_coord)

        if d_of_illumina_sites:
            if f'{chrom}*{orientation}' in d_of_illumina_sites.keys():
                l_correct_pos = find_closest_sites(pos_of_introns, d_of_illumina_sites[f'{chrom}*{orientation}'])
            else:
                l_correct_pos = []
        else:
            l_correct_pos = pos_of_introns

        str_introns = ';'.join(['%s-%s' % pos for pos in l_correct_pos])
        res_line = f'{read_id}\t{chrom}\t{orientation}\t{start_coord}\t{end_coord}\t{str_introns}'
        l_res_lines.append(res_line)
    return l_res_lines


def make_final_file(out_filename, merged_results):
    with open(out_filename, 'w') as ouf:
        ouf.write('\n'.join(merged_results) + '\n')


def change_to_correct_splice_sites(sam_file, ill_sites, num_threads, common_outdir):
    logging.info('    Identify splice sites')
    changed_splice_sites_dir = common_outdir + 'changed_splice_sites/'
    if not os.path.exists(changed_splice_sites_dir):
        os.mkdir(changed_splice_sites_dir)
    out_fname = changed_splice_sites_dir + os.path.basename(sam_file).replace('.sam', '.tsv')

    sam_lines = read_sam_file(sam_file)
    data_by_proc = prepare_data_for_multiprocessing(sam_lines, num_threads)
    if ill_sites:
        d_of_ill_sites = get_d_illumina_splice_sites(ill_sites)
    else:
        d_of_ill_sites = None

    results = Parallel(n_jobs=num_threads)(delayed(extract_splice_sites)(lines_l, d_of_ill_sites)
                                           for lines_l in data_by_proc.values())
    merged_results = functools.reduce(lambda x, y: x.extend(y) or x, results)
    make_final_file(out_fname, merged_results)
    return changed_splice_sites_dir
