import os

from __init__ import PATHS, quiet_call


def align(fastq_1, ref_fn, output_fn, fastq_2=None, p=1):
    '''
    Align reads using Bowtie2.
    '''
    assert os.path.exists(fastq_1)
    if fastq_2:
        assert os.path.exists(fastq_2)

    if fastq_2:
        quiet_call([
            PATHS['bowtie2'],
            '-q',
            '--phred33',
            '-p', str(p),
            '-I', '0',
            '-X', '1000',
            '--fr',
            '--local',
            '--sensitive-local',
            '-S', output_fn,
            '-x', ref_fn,
            '-1', fastq_1,
            '-2', fastq_2,
        ])

    else:
        quiet_call([
            PATHS['bowtie2'],
            '-q',
            '--phred33',
            '-p', str(p),
            '-I', '0',
            '-X', '1000',
            '--local',
            '--sensitive-local',
            '-S', output_fn,
            '-x', ref_fn,
            '-U', fastq_1,
        ])


def build(ref_fn):
    '''
    Build Bowtie2 index for reference if none found.
    '''
    if not (
        os.path.isfile(ref_fn + '.1.bt2') and
        os.path.isfile(ref_fn + '.2.bt2') and
        os.path.isfile(ref_fn + '.3.bt2') and
        os.path.isfile(ref_fn + '.4.bt2') and
        os.path.isfile(ref_fn + '.rev.1.bt2') and
        os.path.isfile(ref_fn + '.rev.2.bt2')
    ):
        quiet_call([
            PATHS['bowtie2_build'],
            ref_fn,
            ref_fn,
        ])
