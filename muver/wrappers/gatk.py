import sys
from subprocess import call

from __init__ import PATHS, quiet_call


#  Run GATK BaseRecalibrator
def run_base_recalibrator(bam, known_sites, ref_fn, recal_table, log_file,):

    quiet_call([
        'java', '-Xmx8g', '-jar',
        PATHS['gatk'],
        '-T', 'BaseRecalibrator',
        '-I', bam,
        '-knownSites', known_sites,
        '--out', recal_table,
        '-R', ref_fn,
        '-cov', 'RepeatLengthCovariate',
        '-cov', 'RepeatUnitCovariate',
        '-log', log_file,
    ])


#  Run GATK PrintReads
def run_print_reads_bqsr(input_bam, ref_fn, recal_table, output_bam, log_file):

    sys.stdout.write('Running GATK PrintReads.\n')
    quiet_call([
        'java', '-Xmx8g', '-jar',
        PATHS['gatk'],
        '-T', 'PrintReads',
        '-I', input_bam,
        '-R', ref_fn,
        '-BQSR', recal_table,
        '-o', output_bam,
        '-log', log_file,
    ])


#  Run GATK HaplotypeCaller
def run_haplotype_caller(bams, ref_fn, output_vcf, log_file, nct=1):

    input_list = []
    for bam in bams:
        input_list.extend(('-I', bam))

    sys.stdout.write('Running GATK HaplotypeCaller.\n')
    quiet_call([
        'java', '-Xmx8g', '-jar',
        PATHS['gatk'],
        '-T', 'HaplotypeCaller',
        '-o', output_vcf,
        '-A', 'StrandAlleleCountsBySample',
        '-A', 'DepthPerSampleHC',
        '-R', ref_fn,
        '-nct', str(nct),
        '-mmq', '5',
        '-log', log_file,
    ] + input_list)


def realigner_target_creator(ref_fn, in_bam, intervals):

    quiet_call([
        'java', '-Xmx4g', '-jar',
        PATHS['gatk'],
        '-R', ref_fn,
        '-T', 'RealignerTargetCreator',
        '-I', in_bam,
        '-o', intervals,
    ])


def indel_realigner(ref_fn, log, in_bam, intervals, realigned_bam):

    quiet_call([
        'java', '-Xmx4g', '-jar',
        PATHS['gatk'],
        '-R', ref_fn,
        '-T', 'IndelRealigner',
        '--maxReadsForRealignment', '100000',
        '-log', log,
        '-I', in_bam,
        '-targetIntervals', intervals,
        '-o', realigned_bam,
    ])
