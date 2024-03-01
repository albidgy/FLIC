import os
from collections import defaultdict


def get_d_starts_or_poly_a(file):
    d_starts_or_poly_a = {}
    with open(file) as starts:
        for line in starts:
            line_l = line.strip('\n').split('\t')
            chrom_and_orientation = f'{line_l[0]}*{line_l[1]}'
            start_coord = int(line_l[2])
            end_coord = int(line_l[3])

            if chrom_and_orientation not in d_starts_or_poly_a.keys():
                d_starts_or_poly_a[chrom_and_orientation] = {}

            for coord in range(start_coord, end_coord + 1):
                d_starts_or_poly_a[chrom_and_orientation][coord] = (start_coord, end_coord)
    return d_starts_or_poly_a


def find_start_stop_for_read(d_cur_chrom_and_orientation, coordinate):
    if coordinate in d_cur_chrom_and_orientation.keys():
        start, stop = d_cur_chrom_and_orientation[coordinate]
        start += 1  # convert 0-based to 1-based
        return f'{start}-{stop}'
    return None


def make_isoforms(splice_sites_file, d_starts, d_poly_a):
    d_of_isoforms = defaultdict(int)
    with open(splice_sites_file) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            chrom = line_l[1]
            orientation = line_l[2]
            start = int(line_l[3])
            stop = int(line_l[4])
            splice_sites = line_l[5]

            if f'{chrom}*{orientation}' not in d_starts.keys():
                continue
            if f'{chrom}*{orientation}' not in d_poly_a.keys():
                continue

            if orientation == '+':
                found_cagefightr_start = find_start_stop_for_read(d_starts[f'{chrom}*{orientation}'], start)
                found_cagefightr_stop = find_start_stop_for_read(d_poly_a[f'{chrom}*{orientation}'], stop)
            elif orientation == '-':
                found_cagefightr_start = find_start_stop_for_read(d_poly_a[f'{chrom}*{orientation}'], start)
                found_cagefightr_stop = find_start_stop_for_read(d_starts[f'{chrom}*{orientation}'], stop)

            if found_cagefightr_start is not None and found_cagefightr_stop is not None:
                isoform = f'{chrom}\t{orientation}\t{found_cagefightr_start}\t{splice_sites}\t{found_cagefightr_stop}'
                d_of_isoforms[isoform] += 1
    return d_of_isoforms


def filter_iso_by_thr(d_of_isoforms, out_filename):
    with open(out_filename, 'w') as ouf:
        for key, val in d_of_isoforms.items():
            ouf.write(f'{key}\t{str(val)}\n')


def create_isoforms(filt_cagefightr_dir, changed_splice_sites_dir, changed_splice_file, common_outdir):
    iso_dir = common_outdir + 'isoforms/'
    if not os.path.exists(iso_dir):
        os.mkdir(iso_dir)

    d_starts = get_d_starts_or_poly_a(f'{filt_cagefightr_dir}filtered_start_cagefightr.bed')
    d_poly_a = get_d_starts_or_poly_a(f'{filt_cagefightr_dir}filtered_polya_cagefightr.bed')
    d_of_isoforms = make_isoforms(changed_splice_sites_dir + changed_splice_file, d_starts, d_poly_a)
    filter_iso_by_thr(d_of_isoforms, f'{iso_dir}isoforms_{changed_splice_file}')
    return iso_dir


def read_isoform_files(path_to_dir, thr_1st, thr_2nd):
    d_of_isoforms = defaultdict(list)
    d_good_iso = {}

    for file in os.listdir(path_to_dir):
        with open(f'{path_to_dir}{file}') as rep:
            for line in rep:
                *line_l, cov = line.strip('\n').split('\t')
                cov = int(cov)
                d_of_isoforms['\t'.join(line_l)].append(cov)

    for key, val_l in d_of_isoforms.items():
        if sum([x >= thr_1st for x in val_l]) > 1 and sum([x >= thr_2nd for x in val_l]) > 0:
            d_good_iso[key] = len(val_l)
    return list(d_good_iso.keys())


def read_genes_file(genes_file):
    d_of_genes = defaultdict(dict)
    with open(genes_file) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            chrom_and_orientation = f'{line_l[0]}*{line_l[1]}'
            start = int(line_l[2])
            stop = int(line_l[4])
            gene_id = line_l[5]

            d_of_genes[chrom_and_orientation][(start, stop)] = gene_id
    return d_of_genes


def intersection_with_genes(iso_coords, d_of_ref):
    start_iso, stop_iso = iso_coords

    for key, gene_id in d_of_ref.items():
        start_ref, stop_ref = key

        max_start = max(start_ref, start_iso)
        min_end = min(stop_ref, stop_iso)
        prop_intersection = max(min_end - max_start, 0) / (stop_iso - start_iso)
        if prop_intersection >= 0.95:
            return gene_id
    return 'None'


def move_isoids_from_genes(iso_dir, genes_file, thr_1st, thr_2nd, out_dir):
    l_of_iso = read_isoform_files(iso_dir, thr_1st, thr_2nd)
    d_of_genes = read_genes_file(genes_file)
    ouf_name = f'{out_dir}final_results/isoforms.tsv'

    d_of_iso_by_genes = defaultdict(list)
    with open(ouf_name, 'w') as ouf:
        ouf.write('')

    for line in l_of_iso:
        line_l = line.split('\t')
        chrom_and_orientation = f'{line_l[0]}*{line_l[1]}'
        start = int(line_l[2].split('-')[0])
        stop = int(line_l[4].split('-')[1])

        gene_id = intersection_with_genes((start, stop), d_of_genes[chrom_and_orientation])
        d_of_iso_by_genes[gene_id].append(((stop - start), line_l))

    for gene, val in d_of_iso_by_genes.items():
        iso_counter = 0
        cur_isoforms = sorted(val)[::-1]
        for elem in cur_isoforms:
            _, iso_l = elem
            iso_counter += 1
            iso_l.append(f'{gene}.{iso_counter}')

            with open(ouf_name, 'a') as ouf:
                ouf.write('\t'.join(iso_l) + '\n')
    return ouf_name
