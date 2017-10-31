import subprocess


def read_excluded_regions(excluded_regions_fn):
    '''
    Read excluded regions from an input file. Return as a dict.
    '''
    excluded_regions = set()

    with open(excluded_regions_fn) as f:
        for line in f:
            chromosome, start, end = line.strip().split()
            start = int(start)
            end = int(end)
            for i in range(start, end + 1):
                excluded_regions.add((
                    chromosome,
                    i,
                ))

    return excluded_regions


def read_cnv_bedgraph(cnv_bedgraph):
        '''
        Read CNV regions from an input bedGraph file. Return as a dict.
        '''
        cnv_regions = dict()

        with open(cnv_bedgraph) as f:

            for line in f:

                chromosome, start, end, ploidy = line.strip().split()
                start = int(start) + 1
                end = int(end)
                ploidy = int(ploidy)

                for i in range(start, end + 1):
                    cnv_regions[(chromosome, i)] = ploidy

        return cnv_regions


def read_chrom_sizes(reference_assembly_fn):
    '''
    Iterate through a reference assembly to find length of associated
    sequences.
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

def read_repeats_var(repeats_fn, var):
    '''
    Read repeats from an input file and return in a dict only those
    that intersect with a given list of variants.
    '''
    repeats = dict()

    with open(repeats_fn) as f:

        for line in f:
            chromosome, sequence, unit_length, unit, start, end = \
                line.strip().split('\t')
            start = int(start)
            end = int(end)
            if chromosome not in repeats:
                repeats[chromosome] = dict()
            if len(sequence) >= 4:
                for i in range(start, end + 1):
                    if (chromosome, i) in var:
                        if i in repeats[chromosome]:
                            repeats[chromosome][i].append({
                                'sequence': sequence,
                                'unit': unit,
                                'start': start,
                            })

                        else:
                            repeats[chromosome][i] = [{
                                'sequence': sequence,
                                'unit': unit,
                                'start': start,
                            }]

    return repeats

def read_repeats(repeats_fn):
    '''
    Read repeats from an input file and return in a dict.
    '''
    repeats = dict()

    with open(repeats_fn) as f:

        for line in f:
            chromosome, sequence, unit_length, unit, start, end = \
                line.strip().split('\t')
            start = int(start)
            end = int(end)
            if chromosome not in repeats:
                repeats[chromosome] = dict()
            if len(sequence) >= 4:
                for i in range(start, end + 1):
                    if i in repeats[chromosome]:
                        repeats[chromosome][i].append({
                            'sequence': sequence,
                            'unit': unit,
                            'start': start,
                        })
                    else:
                        repeats[chromosome][i] = [{
                            'sequence': sequence,
                            'unit': unit,
                            'start': start,
                        }]

    return repeats


def get_mpileup_output(in_bam, ref_fn, out_txt):
    '''
    Pipe the output of samtools mpileup to a TXT file.
    '''
    samtools_path = 'samtools'

    with open(out_txt, 'w') as OUT:
        proc = subprocess.Popen([
            samtools_path, 'mpileup',
            '-q', '5',
            '-Q', '10',
            '-B',
            '-d', '100000',
            '-f', ref_fn,
            in_bam,
        ], stdout=OUT)

        proc.wait()


def read_filtered_sites(samples):
    '''
    Considering a list of sample objects, open the associated filtered sites
    files and add those sites to sample-specific sets. Return the sets in a
    dict.
    '''
    filtered_sites = dict()

    for sample in samples:
        filtered_sites[sample] = set()

        with open(sample.filtered_sites) as f:
            for line in f:
                chromosome, start, end = line.strip().split('\t')
                for position in range(int(start) + 1, int(end) + 1):
                    filtered_sites[sample].add((chromosome, position))

    return filtered_sites
