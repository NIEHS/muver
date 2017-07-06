import math
import numpy
import subprocess


def read_cnv_files(file_list):

    cnv_regions = dict()

    with open(file_list) as _list:
        for entry in _list:

            sample, fn = entry.strip().split()
            ploidy_locus_table[sample] = dict()

            with open(fn) as f:
                for line in f:
                    chromosome, start, end, ploidy = line.strip().split()

                    start = int(start) + 1
                    end = int(end)

                    for i in range(start, end + 1):
                        ploidy_locus_table[sample][(chromosome, str(i))] = \
                            ploidy

    return ploidy_locus_table


def read_excluded_regions(excluded_regions_fn):
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


def read_chrom_sizes(reference_assembly_fn):

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


def read_repeats(repeats_fn):

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


def read_fits(fits_fn):

    fits = dict()

    with open(fits_fn) as f:
        for line in f:
            sample, depth_mean, depth_stdev, strand_bias_stdev = line.strip().split()
            fits[sample] = {
                'depth_mean': float(depth_mean),
                'depth_stdev': float(depth_stdev),
                'strand_bias_stdev': float(strand_bias_stdev),
            }

    return fits


def read_from_distribution(distribution_fn):

    with open(distribution_fn) as f:

        mean = float(f.next().split(':')[1].strip())
        std = float(f.next().split(':')[1].strip())

    return mean, std


def get_mpileup_output(in_bam, ref_fn, out_txt):

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

    filtered_sites = dict()

    for sample in samples:
        filtered_sites[sample] = set()

        with open(sample.filtered_sites) as f:
            for line in f:
                filtered_sites[sample].add(tuple(line.strip().split('\t')[:2]))

    return filtered_sites
