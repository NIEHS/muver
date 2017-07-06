import os
import sys
from subprocess import call

from __init__ import PATHS


#  Align sequences from FASTQ files using Bowtie2
def align(fastq_1, ref_fn, output_fn, fastq_2=None, p=1):

    assert os.path.exists(fastq_1)
    if fastq_2:
        assert os.path.exists(fastq_2)

    if fastq_2:
        call([
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
        call([
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

    if not (
        os.path.isfile(ref_fn + '.1.bt2') and
        os.path.isfile(ref_fn + '.2.bt2') and
        os.path.isfile(ref_fn + '.3.bt2') and
        os.path.isfile(ref_fn + '.4.bt2') and
        os.path.isfile(ref_fn + '.rev.1.bt2') and
        os.path.isfile(ref_fn + '.rev.2.bt2')
    ):
        sys.stdout.write('Bowtie2 index not found: Building.\n')
        call([
            PATHS['bowtie2_build'],
            ref_fn,
            ref_fn,
        ])
