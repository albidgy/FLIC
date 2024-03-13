import logging
from collections import defaultdict


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
            orientation = line_l[1]
            if orientation == '+':
                start = int(line_l[2].split('-')[0])
                end = int(line_l[4].split('-')[1])
            else:
                start = int(line_l[4].split('-')[0])
                end = int(line_l[2].split('-')[1])

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

    with open(f'{out_dir}annotation.gtf', 'w') as ouf:
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
                if orientation == '+':
                    exons_l = get_exons_l(splice_sites, transcript_start, transcript_end)
                else:
                    exons_l = get_exons_l(splice_sites, transcript_end, transcript_start)
                for cur_exon in exons_l:
                    res_line = fill_annot_template(chrom, 'exon', cur_exon[0],
                                                   cur_exon[1], orientation, gene_id,
                                                   transcript_id, exon_number=exon_counter)
                    exon_counter += 1
                    ouf.write(res_line)
