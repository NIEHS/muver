import os

from __init__ import PATHS, quiet_call
from samtools import index_bam


def add_read_groups(in_sam, out_bam, sample_header, tmp_dir, max_records):
    '''
    Run Picard AddOrReplaceReadGroups.

    sample_header -- set read groups to sample_header
    '''
    quiet_call([
        'java', '-jar',
        PATHS['picard'],
        'AddOrReplaceReadGroups',
        'VALIDATION_STRINGENCY=SILENT',
        'SO=coordinate',
        'RGPL=illumina',
        'RGPU=' + sample_header,
        'RGSM=' + sample_header,
        'RGLB=' + sample_header,
        'RGID=' + sample_header,
        'I=' + in_sam,
        'O=' + out_bam,
        'TMP_DIR=' + tmp_dir,
        'MAX_RECORDS_IN_RAM=' + str(max_records),
    ])


def deduplicate(in_bam, out_bam, metrics_file, tmp_dir, max_records):
    '''
    Run Picard MarkDuplicates.
    '''
    quiet_call([
        'java', '-jar',
        PATHS['picard'],
        'MarkDuplicates',
        'VALIDATION_STRINGENCY=SILENT',
        'REMOVE_DUPLICATES=TRUE',
        'I=' + in_bam,
        'O=' + out_bam,
        'M=' + metrics_file,
        'TMP_DIR=' + tmp_dir,
        'MAX_RECORDS_IN_RAM=' + str(max_records),
    ])
    index_bam(out_bam)


def create_sequence_dictionary(ref_fn):
    '''
    Run Picard CreateSequenceDictionary if none found.
    '''
    if not os.path.exists(ref_fn.split('.fa')[0] + '.dict'):
        quiet_call([
            'java', '-jar',
            PATHS['picard'],
            'CreateSequenceDictionary',
            'R=' + ref_fn,
            'O=' + ref_fn.split('.fa')[0] + '.dict',
        ])


def fix_mate_information(in_bam, out_bam, tmp_dir, max_records):
    '''
    Run Picard FixMateInformation.
    '''
    quiet_call([
        'java', '-jar',
        PATHS['picard'],
        'FixMateInformation',
        'VALIDATION_STRINGENCY=SILENT',
        'SO=coordinate',
        'I=' + in_bam,
        'O=' + out_bam,
        'TMP_DIR=' + tmp_dir,
        'MAX_RECORDS_IN_RAM=' + str(max_records),
    ])
