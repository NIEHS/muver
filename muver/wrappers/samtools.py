import os
import sys
import subprocess

from subprocess import call

from __init__ import PATHS, quiet_call


def mapq_filter(in_sam, out_sam, q=20):

    quiet_call([
        PATHS['samtools'], 'view',
        '-Sh',
        '-q', str(q),
        in_sam,
    ], stdout=out_sam)


def index_bam(bam_fn):

    quiet_call([
        PATHS['samtools'], 'index',
        bam_fn,
    ])


def merge_bams(bam_list, out_bam):

    quiet_call([
        PATHS['samtools'], 'merge', '-f',
        out_bam,
    ] + bam_list)

    index_bam(out_bam)


def run_mpileup(bam_file, ref_fn, output_file):
    sys.stdout.write('Running samtools mpileup.\n')
    quiet_call([
        PATHS['samtools'], 'mpileup',
        '-q', '5',
        '-Q', '10',
        '-B',
        '-d', '100000',
        '-f', ref_fn,
        bam_file,
    ], stdout=output_file)


def mpileup_iter(bam_file, ref_fn):

    proc = subprocess.Popen([
        PATHS['samtools'], 'mpileup',
        '-q', '5',
        '-Q', '10',
        '-B',
        '-d', '100000',
        '-f', ref_fn,
        bam_file,
    ], stdout=subprocess.PIPE, stderr=open(os.devnull, 'w'))

    return iter(proc.stdout.readline, '')


def faidx_index(ref_fn):

    if not os.path.exists(ref_fn + '.fai'):
        quiet_call([
            PATHS['samtools'], 'faidx',
            ref_fn,
        ])


def view_bam(input_bam):

    proc = subprocess.Popen([
        PATHS['samtools'], 'view',
        input_bam,
    ], stdout=subprocess.PIPE, stderr=open(os.devnull, 'w'))

    return iter(proc.stdout.readline, '')


def get_mpileup_depths(input_bam, ref_fn, output_bedgraph):

    def print_line(chromosome, start, end, value, OUT):
        if value:  # Only prints non-zero values
            OUT.write(
                '{}\t{}\t{}\t{}\n'.format(
                    chromosome,
                    str(start - 1),
                    str(end),
                    str(value),
                )
            )

    #  Run samtools mpileup
    proc = subprocess.Popen([
        PATHS['samtools'], 'mpileup',
        '-q', '5',
        '-Q', '10',
        '-B',
        '-d', '100000',
        '-f', ref_fn,
        input_bam,
    ], stdout=subprocess.PIPE, stderr=open(os.devnull, 'w'))

    last_pos = None
    last_val = None
    last_chr = None
    start = None
    end = None

    with open(output_bedgraph, 'w') as OUT:
        for line in iter(proc.stdout.readline, ''):
            line_split = line.strip().split()

            chromosome, position, reference_base, coverage = line_split[:4]
            position = int(position)
            coverage = int(coverage)
            if int(coverage) > 0:
                bases = line_split[4]
            else:
                bases = ''

            i = 0
            while i < len(bases):
                if bases[i] == '^':
                    i += 1
                elif bases[i] == '*':
                    coverage += 1
                i += 1

            if chromosome == last_chr and coverage == last_val and position == last_pos + 1:
                end = position
            else:
                print_line(last_chr, start, end, last_val, OUT)
                start = position
                end = position

            last_pos = position
            last_val = coverage
            last_chr = chromosome

        print_line(last_chr, start, end, last_val, OUT)
