import os

from wrappers import bowtie2, picard, samtools


def create_reference_indices(ref_fn):
    '''
    For a given reference FASTA file, generate several indices.
    '''
    bowtie2.build(ref_fn)
    samtools.faidx_index(ref_fn)
    picard.create_sequence_dictionary(ref_fn)


def check_reference_indices(ref_fn):
    '''
    For a given reference FASTA file, check for reference indices.
    '''
    # Bowtie2 and samtools suffixes
    suffixes = [
        '.1.bt2',
        '.2.bt2',
        '.3.bt2',
        '.4.bt2',
        '.rev.1.bt2',
        '.rev.2.bt2',
        '.fai',
    ]
    check_suffixes = \
        [os.path.exists('{}{}'.format(ref_fn, suffix)) for suffix in suffixes]

    picard_dict = '{}{}'.format(os.path.splitext(ref_fn)[0], '.dict')
    check_suffixes.append(os.path.exists(picard_dict))

    return all(check_suffixes)


def read_chrom_sizes(reference_assembly_fn):
    '''
    Iterate through a FASTA file to find the length of each chromosome. If a
    FAIDX index is available, it will read the lengths from there.
    '''
    chrom_sizes = dict()

    if os.path.exists(reference_assembly_fn + '.fai'):
        with open(reference_assembly_fn + '.fai') as f:
            for line in f:
                chromosome, size = line.strip().split('\t')[:2]
                chrom_sizes[chromosome] = int(size)
    else:
        last_chromosome = None
        with open(reference_assembly_fn) as f:
            for line in f:
                if line.startswith('>'):
                    last_chromosome = line.split('>')[1].strip()
                    chrom_sizes[last_chromosome] = 0
                else:
                    chrom_sizes[last_chromosome] += len(line.strip())

    return chrom_sizes
