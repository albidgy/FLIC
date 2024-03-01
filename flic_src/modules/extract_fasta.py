import os
from collections import defaultdict


D_OF_COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '-': '-'}


def complement(seq_l):
    complement_seq_l = []
    for nucl in seq_l:
        if nucl in D_OF_COMPLEMENT.keys():
            complement_seq_l.append(D_OF_COMPLEMENT[nucl])
        else:
            complement_seq_l.append('N')
    return complement_seq_l[::-1]


def read_fasta(fasta_file):
    d_of_fasta = {}
    with open(fasta_file) as fasta:
        for line in fasta:
            if line.startswith('>'):
                chrom = line.split(' ')[0][1:]
                d_of_fasta[chrom] = []
            else:
                line_l = list(line.strip('\n').upper())
                d_of_fasta[chrom].extend(line_l)
    return d_of_fasta


def split_fasta_by_n_nucl(seq):
    seq_splitted = []
    prev_idx = 0
    for idx in range(60, len(seq), 60):
        seq_splitted.append(seq[prev_idx:idx])
        prev_idx = idx
    seq_splitted.append(seq[prev_idx:len(seq)])
    return '\n'.join(seq_splitted)


def change_introns_to_gap(seq_l, splice_sites, start_coord):
    for site in splice_sites:
        start, stop = site
        start_corrected = start - start_coord
        stop_corrected = stop - start_coord

        for idx in range(start_corrected, stop_corrected):
            seq_l[idx] = '-'
    return seq_l


def extract_gene_seqs(file_gene_coords, d_of_fasta, out_dir):
    d_of_gene_seqs = {}
    d_of_gene_coords = {}
    out_fname = out_dir + os.path.basename(file_gene_coords).replace('.tsv', '_seq.fasta')

    with open(file_gene_coords) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            chrom = line_l[0]
            orientation = line_l[1]
            start = int(line_l[2]) - 1
            stop = int(line_l[4])
            gene_id = line_l[5]

            if orientation == '+':
                gene_seq = d_of_fasta[chrom][start:stop]
            elif orientation == '-':
                gene_seq = complement(d_of_fasta[chrom][start:stop])

            d_of_gene_seqs[gene_id] = gene_seq
            d_of_gene_coords[gene_id] = (start, stop)

    with open(out_fname, 'w') as ouf:
        for key, val in d_of_gene_seqs.items():
            ouf.write(f'>{key}\n')
            ouf.write(split_fasta_by_n_nucl(''.join(val)) + '\n')
    return d_of_gene_seqs, d_of_gene_coords


def compare_starts_ends(transcript_coord, gene_coord):
    return abs(gene_coord - transcript_coord)


def extract_transcripts_seq(file_isoform_coords, d_of_gene_coords, d_of_fasta, out_dir):
    d_of_transcripts = defaultdict(dict)
    out_fname = out_dir + os.path.basename(file_isoform_coords).replace('.tsv', '_seq.fasta')
    with open(out_fname, 'w') as ouf:
        ouf.write('')

    with open(file_isoform_coords) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            chrom = line_l[0]
            orientation = line_l[1]
            start = int(line_l[2].split('-')[0]) - 1
            stop = int(line_l[4].split('-')[1])
            transcript_id = line_l[5]
            gene_id = '.'.join(transcript_id.split('.')[:-1])
            d_of_transcripts[gene_id][transcript_id] = []

            if line_l[3] == '':
                splice_sites = line_l[3]
            else:
                splice_sites = [tuple(map(int, x.split('-'))) for x in line_l[3].split(';')]
                splice_sites = sorted([(x[0] - 1, x[1]) for x in splice_sites])  # make 0-based

            transcript_seq = d_of_fasta[chrom][start:stop]
            correct_transcript_seq = change_introns_to_gap(transcript_seq, splice_sites, start)
            # add 5'- and 3'-gaps
            correct_transcript_seq = (['-'] * compare_starts_ends(start, d_of_gene_coords[gene_id][0])) + correct_transcript_seq
            correct_transcript_seq = correct_transcript_seq + (['-'] * compare_starts_ends(stop, d_of_gene_coords[gene_id][1]))

            if orientation == '-':
                correct_transcript_seq = complement(correct_transcript_seq)

            d_of_transcripts[gene_id][transcript_id].extend(correct_transcript_seq)

            with open(out_fname, 'a') as ouf:
                final_seq = split_fasta_by_n_nucl(''.join(correct_transcript_seq).replace('-', ''))
                ouf.write(f'>{transcript_id}\n')
                ouf.write(final_seq + '\n')
    return d_of_transcripts


def make_pseudoalignments(d_of_gene_seqs, d_of_transcripts, out_dir):
    path_dir = out_dir + 'pseudoalignments/'
    os.mkdir(path_dir)

    for gene, gene_seq in d_of_gene_seqs.items():
        with open(f'{path_dir}{gene}.fasta', 'w') as ouf:
            ouf.write(f'>{gene}\n')
            ouf.write(''.join(gene_seq) + '\n')

        for transcript, transcript_seq in d_of_transcripts[gene].items():
            with open(f'{path_dir}{gene}.fasta', 'a') as ouf:
                ouf.write(f'>{transcript}\n')
                ouf.write(''.join(transcript_seq) + '\n')


def extract_fasta(fasta_file, final_fpath_genes, final_fpath_iso, out_dir):
    d_of_fasta = read_fasta(fasta_file)
    d_of_gene_seqs, d_of_gene_coords = extract_gene_seqs(final_fpath_genes, d_of_fasta, out_dir)
    d_of_transcripts = extract_transcripts_seq(final_fpath_iso, d_of_gene_coords, d_of_fasta, out_dir)
    make_pseudoalignments(d_of_gene_seqs, d_of_transcripts, out_dir)
