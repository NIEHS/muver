from wrappers import bowtie2, picard, samtools


def create_reference_indices(ref_fn):
    '''
    For a given refere FASTA file, generate several indices.
    '''
    bowtie2.build(ref_fn)
    samtools.faidx_index(ref_fn)
    picard.create_sequence_dictionary(ref_fn)


def read_chrom_sizes(reference_assembly_fn):
    '''
    Iterate through a FASTA file to find the length of each chromosome.
    '''
    chrom_sizes = dict()
    last_chromosome = None

    with open(reference_assembly_fn) as f:

        for line in f:

            if line.startswith('>'):

                last_chromosome = line.split('>')[1].strip()
                chrom_sizes[last_chromosome] = 0

            else:

                chrom_sizes[last_chromosome] += len(line.strip())

    return chrom_sizes


def read_chrom_sizes_from_file(chrom_sizes_fn):
    '''
    Read chromosome sizes from a UCSC chrom_sizes file.
    '''
    chrom_sizes = dict()

    with open(chrom_sizes_fn) as f:
        for line in f:
            chromosome, size = line.strip().split()
            chrom_sizes[chromosome] = int(size)

    return chrom_sizes
