import logging
import os
import re

from collections import defaultdict
from flic_src.scripts import external_tool_runner


ORF_LEN_RE = re.compile('ORF_len=(\d+)')
TRANSCRIPT_ID_RE = re.compile(r'transcript_id "(.+?)"')


def read_genes_file(genes_fpath):
    d_genes = {}
    with open(genes_fpath) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            chrom = line_l[0]
            orientation = line_l[1]
            start = int(line_l[2])
            end = int(line_l[4])
            gene_id = line_l[-1]
            d_genes[gene_id] = (chrom, start, end, orientation)

    return dict(sorted(d_genes.items(), key=lambda x: (x[1], x[0])))  # sort by str(chrom)!


def read_iso_file(isoform_fpath):
    d_iso_comb_by_genes = defaultdict(list)
    with open(isoform_fpath) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            start = int(line_l[2].split('-')[0])
            end = int(line_l[4].split('-')[1])

            splice_sites = line_l[3]
            transcript_id = line_l[-1]
            gene_id = '.'.join(transcript_id.split('.')[:-1])

            d_iso_comb_by_genes[gene_id].append((transcript_id, start, end, splice_sites))

    return d_iso_comb_by_genes


def prepare_introns_l(introns_str):
    correct_introns_l = []
    if introns_str == '':
        return correct_introns_l
    introns_l = introns_str.split(';')

    for elem in introns_l:
        elem = tuple(map(int, elem.split('-')))
        correct_introns_l.append(elem)

    return correct_introns_l


def get_exons_l(introns_str, start_transcript, end_transcript):
    introns_l = prepare_introns_l(introns_str)
    exons_l = []
    next_exon_start = start_transcript

    for intron in introns_l:
        exons_l.append((next_exon_start, intron[0] - 1))
        next_exon_start = intron[1] + 1
    exons_l.append((next_exon_start, end_transcript))

    return exons_l


def fill_annot_template(chrom, feat_type, start, end, orientation,
                        gene_id, transcript_id='', tool='FLIC', exon_number=None):
    if feat_type == 'exon':
        comment_line = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "{str(exon_number)}";'
    else:
        comment_line = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'

    return f'{chrom}\t{tool}\t{feat_type}\t{str(start)}\t{str(end)}\t.\t{orientation}\t.\t{comment_line}\n'


def create_gtf_file(isoform_fpath, genes_fpath, out_dir, version):
    logging.info('Create annotation GTF format file')
    d_genes = read_genes_file(genes_fpath)
    d_iso_comb_by_genes = read_iso_file(isoform_fpath)
    anno_ouf_path = os.path.join(out_dir, 'annotation_wo_cds.gtf')

    with open(anno_ouf_path, 'w') as ouf:
        ouf.write(f'#FLIC {version} created GTF\n')
        for gene_id in d_genes.keys():
            chrom, gene_start, gene_end, orientation = d_genes[gene_id]
            res_line = fill_annot_template(chrom, 'gene', gene_start,
                                           gene_end, orientation, gene_id)
            ouf.write(res_line)

            for transcript_info in d_iso_comb_by_genes[gene_id]:
                transcript_id, transcript_start, transcript_end, splice_sites = transcript_info
                res_line = fill_annot_template(chrom, 'transcript', transcript_start, transcript_end,
                                               orientation, gene_id, transcript_id)
                ouf.write(res_line)

                exon_counter = 1
                exons_l = get_exons_l(splice_sites, transcript_start, transcript_end)
                for cur_exon in exons_l:
                    res_line = fill_annot_template(chrom, 'exon', cur_exon[0],
                                                   cur_exon[1], orientation, gene_id,
                                                   transcript_id, exon_number=exon_counter)
                    exon_counter += 1
                    ouf.write(res_line)

    return anno_ouf_path


def run_orfipy(iso_seq_fpath, num_threads, common_outdir):
    logging.info('    Predict ORFs using orfipy')
    orfipy_outdir = os.path.join(common_outdir, 'predicted_orfs/')
    cmd_orfipy = f'orfipy --min 300 --max 100000 --include-stop --strand f --procs {num_threads} --start ATG --outdir {orfipy_outdir} --bed isoforms_orfs.bed {iso_seq_fpath}'
    external_tool_runner.run_external_tool(cmd_orfipy, common_outdir)
    return os.path.join(orfipy_outdir, 'isoforms_orfs.bed')


def keep_longest_orfs(orfipy_fpath):
    d_longest_orfs_by_iso = {}
    orfipy_oufpath = orfipy_fpath.replace('.bed', '_longest.bed')
    with open(orfipy_fpath) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            iso_id = line_l[0]
            orf_len = int(ORF_LEN_RE.search(line_l[3]).group(1))
            if iso_id not in d_longest_orfs_by_iso:
                d_longest_orfs_by_iso[iso_id] = (0, None)

            if orf_len > d_longest_orfs_by_iso[iso_id][0]:
                d_longest_orfs_by_iso[iso_id] = (orf_len, line)

    with open(orfipy_oufpath, 'w') as ouf:
        for iso_id in d_longest_orfs_by_iso:
            ouf.write(d_longest_orfs_by_iso[iso_id][1])

    return orfipy_oufpath


def read_bed_file(bed_fpath):
    d_orf_coords = {}

    with open(bed_fpath) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            iso_id = line_l[0]
            orf_start = int(line_l[1]) + 1  # convert to 1-based
            orf_end = int(line_l[2])
            d_orf_coords[iso_id] = (orf_start, orf_end)

    return d_orf_coords


def extract_exon_coords(anno_fpath):
    d_iso_info = {}

    with open(anno_fpath) as gtf:
        for line in gtf:
            if line[0] == '#':
                continue
            line_l = line.strip('\n').split('\t')

            if line_l[2] == 'exon':
                iso_id = TRANSCRIPT_ID_RE.search(line_l[8]).group(1)
                start = int(line_l[3])
                end = int(line_l[4])
                strand = line_l[6]
                if iso_id not in d_iso_info:
                    d_iso_info[iso_id] = []
                d_iso_info[iso_id].append((start, end, strand))

    for iso_id in d_iso_info:
        d_iso_info[iso_id] = sorted(d_iso_info[iso_id])

    return d_iso_info


def make_corr_cds_coords(orf_coords, iso_struct):
    cds_coords = []
    orf_start, orf_end = orf_coords
    strand = iso_struct[0][2]
    transcript_len = sum(end - start + 1 for start, end, _ in iso_struct)
    transcript_pos = 1

    if strand == '-':
        orf_start, orf_end = (transcript_len - orf_end + 1, transcript_len - orf_start + 1)

    for exon_start, exon_end, _ in iso_struct:
        exon_len = exon_end - exon_start + 1
        exon_t_start = transcript_pos
        exon_t_end = transcript_pos + exon_len - 1

        overlap_start = max(orf_start, exon_t_start)
        overlap_end = min(orf_end, exon_t_end)
        if overlap_start <= overlap_end:
            g_start = exon_start + (overlap_start - exon_t_start)
            g_end = exon_start + (overlap_end - exon_t_start)
            cds_coords.append([g_start, g_end])
        transcript_pos += exon_len
        if transcript_pos > orf_end:
            break

    return cds_coords


def add_cds_lines(anno_fpath, d_cds_coords, anno_oufpath):
    d_anno_prev = {}
    d_cds_supp_info = {}
    with open(anno_fpath) as inf, open(anno_oufpath, 'w') as ouf:
        for line in inf:
            if line[0] == '#':
                ouf.write(line)
                continue
            line_l = line.strip('\n').split('\t')
            if line_l[2] == 'gene':
                gene_line = line
                is_recorded_gene = False
            else:
                transcript_id = TRANSCRIPT_ID_RE.search(line_l[8]).group(1)
                if transcript_id not in d_anno_prev:
                    d_anno_prev[transcript_id] = []
                    d_cds_supp_info[transcript_id] = [line_l[0], line_l[1], line_l[6], line_l[8]]

                if not is_recorded_gene:
                    d_anno_prev[transcript_id].append(gene_line)
                    is_recorded_gene = True
                d_anno_prev[transcript_id].append(line)

        for transcript_id in d_anno_prev:
            ouf.write(''.join(d_anno_prev[transcript_id]))
            if transcript_id in d_cds_coords:  # some isoforms may be non-protein coding
                coding_len = 0
                sup_info = d_cds_supp_info[transcript_id]
                for cds_border in d_cds_coords[transcript_id]:
                    phase = (3 - coding_len % 3) % 3
                    coding_len += cds_border[1] - cds_border[0] + 1
                    ouf.write(f'{sup_info[0]}\t{sup_info[1]}\tCDS\t{cds_border[0]}\t{cds_border[1]}\t.\t{sup_info[2]}\t{phase}\t{sup_info[3]}\n')


def add_cds_into_anno(anno_fpath, iso_seq_fpath, num_threads, common_outdir):
    anno_oufpath = os.path.join(common_outdir, os.path.basename(anno_fpath).replace('_wo_cds.gtf', '.gtf'))
    orfipy_fpath = run_orfipy(iso_seq_fpath, num_threads, common_outdir)
    bed_fpath = keep_longest_orfs(orfipy_fpath)
    d_orf_coords = read_bed_file(bed_fpath)
    d_iso_info = extract_exon_coords(anno_fpath)

    d_cds_coords = {}
    for iso_id in d_orf_coords:
        cds_coords = make_corr_cds_coords(d_orf_coords[iso_id], d_iso_info[iso_id])
        d_cds_coords[iso_id] = cds_coords
    add_cds_lines(anno_fpath, d_cds_coords, anno_oufpath)

    os.remove(anno_fpath)  # remove GTF wo CDS
    return os.path.dirname(bed_fpath)
