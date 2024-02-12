import os
import shutil
import sys

from scripts import option_parser
from scripts.external_tool_runner import Logger
from utils import trim_map_ont, downsampling_by_annot, make_correct_splice_sites, find_start_polya, make_iso, \
    make_genes_from_iso, extract_fasta


def run_pipeline():
    arguments = option_parser.parser_arguments()
    if os.path.exists(arguments.output_dir):
        shutil.rmtree(arguments.output_dir)
    os.mkdir(arguments.output_dir)
    if arguments.output_dir[-1] != '/':
        arguments.output_dir = arguments.output_dir + '/'
    sys.stdout = Logger(arguments.output_dir)

    l_of_ont_reads = arguments.nanopore_reads.split(',')

    for long_reads in l_of_ont_reads:
        cur_file_location = os.path.dirname(long_reads) + '/'
        long_reads = os.path.basename(long_reads)

        if arguments.trim_long_reads:
            cur_file_location = trim_map_ont.run_porechop(f'{cur_file_location}{long_reads}',
                                                          arguments.output_dir, arguments.threads)
        if arguments.cut_adapt:
            cur_file_location = trim_map_ont.run_cutadapt(f'{cur_file_location}{long_reads}',
                                                          arguments.output_dir, arguments.threads)

        cur_file_location = trim_map_ont.run_minimap2(f'{cur_file_location}{long_reads}', arguments.ref_fasta,
                                                      arguments.output_dir, arguments.threads)

        long_reads = long_reads.split('.fastq')[0] + '.sam'
        if arguments.make_downsampling:
            cur_file_location, uniq_map_sam = downsampling_by_annot.run_downsampling(f'{cur_file_location}{long_reads}',
                                                                       arguments.ref_annot,
                                                                       arguments.threads,
                                                                       arguments.downsampling_max_thr,
                                                                       arguments.downsampling_min_thr,
                                                                       arguments.output_dir)

        changed_splice_sites_dir = make_correct_splice_sites.change_to_correct_splice_sites(uniq_map_sam,
                                                                                            arguments.ill_sites,
                                                                                            arguments.threads,
                                                                                            arguments.output_dir)
        bam_dir = trim_map_ont.convert_sam2bam(f"{cur_file_location}{long_reads}",
                                               arguments.output_dir, arguments.threads)

    cagefightr_res_dir = find_start_polya.find_starts_polya(bam_dir, arguments.ref_fasta,
                                                            arguments.output_dir)
    # cagefightr_res_dir = '/mnt/nvme/a_kasianova/new_splicing_pipeline/meristem_sampling/0_run_16_mln/peak_calling/final_cagefightr_res/'

    for changed_splice_file in os.listdir(changed_splice_sites_dir):
        iso_dir = make_iso.create_isoforms(cagefightr_res_dir, changed_splice_sites_dir,
                                           changed_splice_file, arguments.output_dir)

    genes_dir = make_genes_from_iso.create_genes(iso_dir, arguments.threads, arguments.thr1,
                                                 arguments.thr2, arguments.output_dir)
    final_fpath_genes = make_genes_from_iso.move_geneids_from_annot(arguments.ref_annot,
                                                                    f'{genes_dir}genes.tsv', arguments.output_dir)
    final_fpath_iso = make_iso.move_isoids_from_genes(iso_dir, final_fpath_genes, arguments.thr1,
                                                      arguments.thr2, arguments.output_dir)
    extract_fasta.extract_fasta(arguments.ref_fasta, final_fpath_genes, final_fpath_iso, arguments.output_dir)


if __name__ == '__main__':
    run_pipeline()
