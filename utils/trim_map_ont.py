import os

from scripts import option_parser, external_tool_runner


def run_porechop(long_reads, common_outdir, num_threads):
    porechop_outdir = common_outdir + 'porechop_output/'
    if not os.path.exists(porechop_outdir):
        os.mkdir(porechop_outdir)
    cmd_porechop = f'{option_parser.fpath_porechop} -t {str(num_threads)} -i {long_reads} -o {porechop_outdir}{os.path.basename(long_reads)}'
    external_tool_runner.run_external_tool(cmd_porechop, common_outdir)
    return porechop_outdir


def run_cutadapt(long_reads, common_outdir, num_threads):
    cutadapt_outdir = common_outdir + 'cutadapt_output/'
    if not os.path.exists(cutadapt_outdir):
        os.mkdir(cutadapt_outdir)
    n_polya_tail = '"A{100}"'
    cmd_cutadapt = f'{option_parser.fpath_cutadapt} -a {n_polya_tail} -j {num_threads} -o {cutadapt_outdir}{os.path.basename(long_reads)} {long_reads}'
    external_tool_runner.run_external_tool(cmd_cutadapt, common_outdir)
    return cutadapt_outdir


def run_minimap2(long_reads, ref_fasta, common_outdir, num_threads):
    minimap2_outdir = common_outdir + 'minimap2_output/'
    if not os.path.exists(minimap2_outdir):
        os.mkdir(minimap2_outdir)
    cmd_minimap2 = f"{option_parser.fpath_minimap2} -ax splice -k14 -uf -t {num_threads} -G 10k {ref_fasta} {long_reads} > {minimap2_outdir}{os.path.basename(long_reads).split('.fastq')[0]}.sam"
    external_tool_runner.run_external_tool(cmd_minimap2, common_outdir)
    return minimap2_outdir


def convert_sam2bam(sam_file, common_outdir, num_threads):
    bam_outdir = common_outdir + 'sorted_bams/'
    if not os.path.exists(bam_outdir):
        os.mkdir(bam_outdir)
    cmd_samtools = f"{option_parser.fpath_samtools} view -@ {num_threads} -bS {sam_file} | samtools sort -@ {num_threads} > {bam_outdir}{os.path.basename(sam_file).replace('.sam', '.sorted.bam')}"
    external_tool_runner.run_external_tool(cmd_samtools, common_outdir)
    return bam_outdir
