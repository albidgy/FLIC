import logging
import os

from flic_src.scripts import external_tool_runner


def run_porechop(long_reads, common_outdir, num_threads):
    logging.info('    Trim long reads using Porechop')
    porechop_outdir = common_outdir + 'porechop_output/'
    if not os.path.exists(porechop_outdir):
        os.mkdir(porechop_outdir)
    cmd_porechop = f'porechop -t {str(num_threads)} -i {long_reads} -o {porechop_outdir}{os.path.basename(long_reads)}'
    external_tool_runner.run_external_tool(cmd_porechop, common_outdir)
    return porechop_outdir


def run_cutadapt(long_reads, common_outdir, num_threads):
    logging.info('    Cut adapters using cutadapt')
    cutadapt_outdir = common_outdir + 'cutadapt_output/'
    if not os.path.exists(cutadapt_outdir):
        os.mkdir(cutadapt_outdir)
    cmd_cutadapt = f'cutadapt -a "A{{100}}" -j {num_threads} -o {cutadapt_outdir}{os.path.basename(long_reads)} {long_reads}'
    external_tool_runner.run_external_tool(cmd_cutadapt, common_outdir)
    return cutadapt_outdir


def run_minimap2(long_reads, ref_fasta, common_outdir, max_intron_len, num_threads):
    logging.info('    Map long reads to reference using minimap2')
    minimap2_outdir = common_outdir + 'minimap2_output/'
    if not os.path.exists(minimap2_outdir):
        os.mkdir(minimap2_outdir)
    cmd_minimap2 = f"minimap2 -ax splice -k14 -uf -t {num_threads} -G {max_intron_len} {ref_fasta} {long_reads} > {minimap2_outdir}{os.path.basename(os.path.splitext(long_reads.split('.gz')[0])[0])}.sam"
    external_tool_runner.run_external_tool(cmd_minimap2, common_outdir)
    return minimap2_outdir


def convert_sam2bam(sam_file, common_outdir, num_threads):
    logging.info('    Convert sam to bam and sort them using samtools\n')
    bam_outdir = common_outdir + 'sorted_bams/'
    if not os.path.exists(bam_outdir):
        os.mkdir(bam_outdir)
    cmd_samtools = f"samtools view -@ {num_threads} -bS {sam_file} | samtools sort -@ {num_threads} > {bam_outdir}{os.path.basename(sam_file).replace('.sam', '.sorted.bam')}"
    external_tool_runner.run_external_tool(cmd_samtools, common_outdir)
    return bam_outdir
