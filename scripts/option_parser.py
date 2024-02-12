import argparse
from datetime import datetime


def parser_arguments():
    parser = argparse.ArgumentParser(description='FLIC: tool for isoform reconstruction based on long reads')
    general_args = parser.add_argument_group(description='General arguments')
    general_args.add_argument('--long_reads',
                              required=True,
                              help='Long reads in fastq or fastq.gz format separated by commas [Example: '
                                   '/path/to/ont_rep1.fastq,/path/to/ont_rep2.fastq]',
                              )
    general_args.add_argument('--ref_fasta',
                              required=True,
                              help='Path to reference FASTA file',
                              )
    general_args.add_argument('--ref_annot',
                              required=True,
                              help='Path to reference annotation file in GTF format',
                              )
    general_args.add_argument('-t',
                              '--threads',
                              default=1,
                              type=int,
                              help='Number of threads [default: 1]',
                              )
    general_args.add_argument('-o',
                              '--output_dir',
                              default='../res_splicing_' + datetime.now().strftime('%Y_%m_%d_%H_%M_%S') + '/',
                              help='Output directory [default: ../res_splicing_YEAR_MONTH_DAY_HOUR_MINUTE_SECOND]',
                              )
    general_args.add_argument('--trim_long_reads',
                              default=False,
                              action='store_true',
                              help='Add trimming long reads step by using Porechop tool [default: False]',
                              )
    general_args.add_argument('--cut_adapt',
                              default=False,
                              action='store_true',
                              help='Cut polyA tail of long reads by using cutadapt tool [default: False]',
                              )
    general_args.add_argument('--make_downsampling',
                              default=False,
                              action='store_true',
                              help='Make a downsampling long reads by given threshold [default: False]',
                              )

    optional_args = parser.add_argument_group(description='Optional arguments')
    optional_args.add_argument('--splice_sites',
                               default=None,
                               help='Path to file with a list of splice sites',
                               )
    optional_args.add_argument('--downsampling_min_thr',
                               default=5,
                               type=int,
                               help='The number of reads per gene less than a specified value is considered noise and '
                                    'is excluded from the analysis [default: 5]',
                               )
    optional_args.add_argument('--downsampling_max_thr',
                               default=1000,
                               type=int,
                               help='Number of reads per gene remaining after downsampling [default: 1000]',
                               )
    optional_args.add_argument('--iso_thr1',  # '--n_reads_per_iso_thr_2nd',
                               default=5,
                               type=int,
                               help='Minimum number of reads forming an isoform at least in 1 replicate [default: 5]'
                               )
    optional_args.add_argument('--iso_thr2',  # '--n_reads_per_iso_thr_1st',
                               default=1,
                               type=int,
                               help='Minimum number of reads forming an isoform in another replicate [default: 1]'
                               )
    return parser.parse_args()


fpath_porechop = 'porechop'
fpath_cutadapt = 'cutadapt'
fpath_minimap2 = 'minimap2'
fpath_samtools = 'samtools'
fpath_bedtools = 'bedtools'
fpath_bedgraph2bigwig = 'bedGraphToBigWig'
