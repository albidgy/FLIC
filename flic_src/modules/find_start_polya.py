import os

from flic_src.scripts import external_tool_runner


D_FOR_BEDTOOLS = {('5', '+'): 'start_fwd', ('5', '-'): 'start_rev',
                  ('3', '+'): 'polya_fwd', ('3', '-'): 'polya_rev'}


def make_fasta_chr_len(ref_fasta, ouf_dir):
    d_of_chr_len = {}
    ouf_name = ouf_dir + 'ref_chr_length.txt'

    with open(ref_fasta) as inf:
        for line in inf:
            if line.startswith('>'):
                chrom = line.strip('\n').split(' ')[0][1:]
                d_of_chr_len[chrom] = 0
            else:
                n_nucls = len(line.strip('\n'))
                d_of_chr_len[chrom] += n_nucls

    with open(ouf_name, 'w') as ouf:
        for key, val in d_of_chr_len.items():
            ouf.write(f'{key}\t{val}\n')
    return ouf_name


def bam2bedgraph(path_to_cutted_bams, peaks_dir):
    output_dirname = peaks_dir + 'bedgraph_files/'
    os.mkdir(output_dirname)

    for file in os.listdir(path_to_cutted_bams):
        out_fname = file.replace('.bam', '.bedgraph')
        for key, val in D_FOR_BEDTOOLS.items():
            cmd_bam2bedgraph = f'bedtools genomecov -ibam {path_to_cutted_bams}{file} -{key[0]} -strand {key[1]} -bg > {output_dirname}{val}_{out_fname}'
            external_tool_runner.run_external_tool(cmd_bam2bedgraph, f"{os.path.split(peaks_dir[:-1])[0]}/")
    return output_dirname


def run_split_bedgraph(file_bedgraph, peaks_dir):
    output_name = peaks_dir + os.path.basename(file_bedgraph)
    l_lines = []

    with open(file_bedgraph) as bedgraph:
        for line in bedgraph:
            chrom, start, stop, cov = line.strip('\n').split('\t')
            start = int(start)
            stop = int(stop)
            while start < stop:
                l_lines.append((chrom, start, start + 1, cov))
                start += 1
    l_lines = sorted(l_lines)
    with open(output_name, 'w') as ouf:
        for elem in l_lines:
            ouf.write('\t'.join(list(map(str, elem))) + '\n')


def make_correct_bedgraph(path_to_bedgraph, peaks_dir):
    output_dirname = peaks_dir + 'corrected_bedgraph/'
    os.mkdir(output_dirname)

    for file in os.listdir(path_to_bedgraph):
        run_split_bedgraph(path_to_bedgraph + file, output_dirname)
    return output_dirname


def bedgraph2bigwig(path_to_corrected_bedgraph, file_chr_length, peaks_dir):
    output_dirname = peaks_dir + 'bigwig_files/'
    os.mkdir(output_dirname)

    for file in os.listdir(path_to_corrected_bedgraph):
        out_fname = file.split('.bedgraph')[0] + '.bw'
        cmd_bedgraph2biwig = f'bedGraphToBigWig {path_to_corrected_bedgraph}{file} {file_chr_length} {output_dirname}{out_fname}'
        external_tool_runner.run_external_tool(cmd_bedgraph2biwig, f"{os.path.split(peaks_dir[:-1])[0]}/")
    return output_dirname


def run_cagefightr(path_to_bigwig, file_chr_length, peaks_dir):
    fpath_cagefightr = os.path.abspath(os.path.dirname(__file__)) + '/run_cagefightr.R'
    d_of_clust_files = {'start_fwd': [], 'start_rev': [],
                        'polya_fwd': [], 'polya_rev': []}
    output_dirname = peaks_dir + 'cagefightr_out/'
    os.mkdir(output_dirname)

    for file in sorted(os.listdir(path_to_bigwig)):
        for key in d_of_clust_files.keys():
            if key in file:
                d_of_clust_files[key].append(path_to_bigwig + file)

    for key, val in d_of_clust_files.items():
        d_of_clust_files[key] = ','.join(val)

    cmd_cagefightr_fwd = f'Rscript {fpath_cagefightr} --forward {d_of_clust_files["start_fwd"]} --reverse {d_of_clust_files["start_rev"]} --chromosome {file_chr_length} --output {output_dirname}start_cagefightr.bed'
    cmd_cagefightr_rev = f'Rscript {fpath_cagefightr} --forward {d_of_clust_files["polya_fwd"]} --reverse {d_of_clust_files["polya_rev"]} --chromosome {file_chr_length} --output {output_dirname}polya_cagefightr.bed'
    external_tool_runner.run_external_tool(cmd_cagefightr_fwd, f"{os.path.split(peaks_dir[:-1])[0]}/")
    external_tool_runner.run_external_tool(cmd_cagefightr_rev, f"{os.path.split(peaks_dir[:-1])[0]}/")
    return output_dirname


def filter_cagefightr(path_to_cagefightr_dir, peaks_dir):
    output_dirname = peaks_dir + 'final_cagefightr_res/'
    os.mkdir(output_dirname)

    for file in os.listdir(path_to_cagefightr_dir):
        with open(f'{output_dirname}filtered_{file}', 'w') as ouf:
            ouf.write('')

        with open(f'{path_to_cagefightr_dir}{file}') as inf:
            for line in inf:
                line_l = line.strip('\n').split('\t')
                chrom = line_l[0]
                start = line_l[1]
                end = line_l[2]
                orientation = line_l[5]

                with open(f'{output_dirname}filtered_{file}', 'a') as ouf:
                    res_line = f'{chrom}\t{orientation}\t{start}\t{end}\n'
                    ouf.write(res_line)
    return output_dirname


def find_starts_polya(bam_dir, ref_fasta, common_outdir):
    peaks_dir = common_outdir + 'peak_calling/'
    if not os.path.exists(peaks_dir):
        os.mkdir(peaks_dir)

    file_chr_length = make_fasta_chr_len(ref_fasta, common_outdir)
    path_to_bedgraph = bam2bedgraph(bam_dir, peaks_dir)
    path_to_corrected_bedgraph = make_correct_bedgraph(path_to_bedgraph, peaks_dir)
    path_to_bigwig = bedgraph2bigwig(path_to_corrected_bedgraph, file_chr_length, peaks_dir)
    path_to_cagefightr_res = run_cagefightr(path_to_bigwig, file_chr_length, peaks_dir)
    path_to_filt_cagefightr_res = filter_cagefightr(path_to_cagefightr_res, peaks_dir)
    os.remove(f'{common_outdir}ref_chr_length.txt')

    return path_to_filt_cagefightr_res
