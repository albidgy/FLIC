import logging
import os
import shutil
import sys

from flic_src.modules import make_iso, trim_map_ont, extract_fasta, downsampling, find_start_polya, \
    make_correct_splice_sites, make_genes_from_iso, filter_by_1percent, create_annot
from flic_src.scripts.external_tool_runner import Logger
from flic_src.scripts import option_parser


def run_tool():
    arguments = option_parser.parser_arguments()
    if os.path.exists(arguments.output_dir):
        shutil.rmtree(arguments.output_dir)
    os.mkdir(arguments.output_dir)
    if arguments.output_dir[-1] != '/':
        arguments.output_dir = arguments.output_dir + '/'

    sys.stdout = Logger(arguments.output_dir)
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s:  %(asctime)s  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S",
                        handlers=[logging.FileHandler(f"{arguments.output_dir}FLIC.log"), logging.StreamHandler()])
    logging.info(f'Command: {"".join(sys.argv)}')
    logging.info('____Getting started FLIC____')

    l_of_ont_reads = arguments.long_reads.split(',')
    if len(l_of_ont_reads) <= 1:
        print('Error: FLIC requires at least 2 repetitions of long reads as input data', file=sys.stderr)
        sys.exit(1)
    if arguments.ref_annot is not None and arguments.ref_annot.split('.')[-1] != 'gtf':
        print('Error: FLIC can work only with annotations in GTF format', file=sys.stderr)
        sys.exit(1)

    dirs_for_delete = set()
    counter_files = 1

    for long_reads in l_of_ont_reads:
        logging.info(f'  File {str(counter_files)}:')
        cur_file_location = os.path.abspath(os.path.dirname(long_reads)) + '/'
        long_reads = os.path.basename(long_reads)

        if long_reads.split('.')[-1] != 'sam':
            if arguments.trim_long_reads:
                cur_file_location = trim_map_ont.run_porechop(f'{cur_file_location}{long_reads}',
                                                              arguments.output_dir, arguments.threads)
                dirs_for_delete.add(cur_file_location)

            if arguments.cut_adapt:
                cur_file_location = trim_map_ont.run_cutadapt(f'{cur_file_location}{long_reads}',
                                                              arguments.output_dir, arguments.threads)
                dirs_for_delete.add(cur_file_location)

            cur_file_location = trim_map_ont.run_minimap2(f'{cur_file_location}{long_reads}', arguments.ref_fasta,
                                                          arguments.output_dir, arguments.threads)
            dirs_for_delete.add(cur_file_location)
            long_reads = os.path.splitext(long_reads.split('.gz')[0])[0] + '.sam'

        cur_file_location, uniq_map_sam = downsampling.run_downsampling(arguments.make_downsampling,
                                                                        arguments.ref_annot,
                                                                        f'{cur_file_location}{long_reads}',
                                                                        arguments.ref_annot,
                                                                        arguments.threads,
                                                                        arguments.downsampling_max_thr,
                                                                        arguments.downsampling_min_thr,
                                                                        arguments.output_dir)
        dirs_for_delete.add(cur_file_location)

        changed_splice_sites_dir = make_correct_splice_sites.change_to_correct_splice_sites(uniq_map_sam,
                                                                                            arguments.splice_sites,
                                                                                            arguments.threads,
                                                                                            arguments.output_dir)
        dirs_for_delete.add(changed_splice_sites_dir)

        if arguments.make_downsampling:
            sample_fname = long_reads
        else:
            sample_fname = 'unique_map_' + long_reads

        if arguments.peaks_dir is None:
            bam_dir = trim_map_ont.convert_sam2bam(f"{cur_file_location}{sample_fname}",
                                                   arguments.output_dir, arguments.threads)
            dirs_for_delete.add(bam_dir)
        counter_files += 1

    if arguments.peaks_dir is None:
        cagefightr_res_dir = find_start_polya.find_starts_polya(bam_dir, arguments.ref_fasta,
                                                                arguments.output_dir)
        dirs_for_delete.add(os.path.split(os.path.abspath(cagefightr_res_dir))[0])
    else:
        cagefightr_res_dir = arguments.peaks_dir

    logging.info('Reconstruct isoforms')
    for changed_splice_file in os.listdir(changed_splice_sites_dir):
        iso_dir = make_iso.create_isoforms(cagefightr_res_dir, changed_splice_sites_dir,
                                           changed_splice_file, arguments.output_dir)
        dirs_for_delete.add(iso_dir)

    genes_dir = make_genes_from_iso.create_genes(iso_dir, arguments.threads, arguments.iso_thr2,
                                                 arguments.iso_thr1, arguments.output_dir)
    dirs_for_delete.add(genes_dir)

    if arguments.ref_annot is not None:
        final_fpath_genes = make_genes_from_iso.move_geneids_from_annot(arguments.ref_annot,
                                                                        f'{genes_dir}genes.tsv',
                                                                        arguments.output_dir)
    else:
        final_fpath_genes = make_genes_from_iso.create_gene_ids_wo_annot(f'{genes_dir}genes.tsv',
                                                                         arguments.output_dir)

    final_fpath_iso = make_iso.move_isoids_from_genes(iso_dir, final_fpath_genes, arguments.iso_thr2,
                                                      arguments.iso_thr1, arguments.output_dir)
    extract_fasta.extract_fasta(arguments.ref_fasta, final_fpath_genes, final_fpath_iso, arguments.output_dir)
    create_annot.create_gtf_file(final_fpath_iso, final_fpath_genes, arguments.output_dir, option_parser.VERSION)

    if arguments.extra_filter_iso:
        filter_by_1percent.filter_by_1percent(iso_dir, arguments.iso_thr1, arguments.iso_thr2,
                                              final_fpath_iso, arguments.output_dir)

    logging.info('Clear intermediate dirs')
    for dir_name in dirs_for_delete:
        shutil.rmtree(dir_name)
    logging.info('____FLIC finished successfully____')


if __name__ == '__main__':
    run_tool()
