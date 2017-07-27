from collections import defaultdict
import re
import numpy
import csv

from wrappers.samtools import mpileup_iter


def get_allelic_fractions(bam_file, ref_fn, output_file):
    '''
    Read over a BAM file, output an allele fraction histogram.
    '''
    allele_fractions = []

    for line in mpileup_iter(bam_file, ref_fn):

        allele_counts = defaultdict(int)

        line_split = line.strip().split()
        reference_allele, coverage = line_split[2:4]
        if int(coverage) > 0:
            bases = line_split[4]
        else:
            bases = ''

        i = 0
        while i < len(bases):
            if re.match('[.,]', bases[i]):
                allele_counts[reference_allele] += 1
            elif re.match('[ACGTacgt]', bases[i]):
                allele_counts[bases[i].upper()] += 1
            elif re.match('[+-]', bases[i]):
                indel_type = bases[i]
                i += 1
                indel_length = int(bases[i])
                i += 1
                indel = indel_type + bases[i:i+indel_length].upper()
                allele_counts[indel] += 1
                i += indel_length
            elif bases[i] == '^':
                i += 1
            elif bases[i] == '*':
                allele_counts[bases[i]] += 1
            i += 1

        _sum = sum(allele_counts.values())
        if _sum > 20:
            for value in allele_counts.values():
                allele_fractions.append(float(value) / _sum)

    bins = [(0.01 * i) for i in range(101)]
    histogram = numpy.histogram(allele_fractions, bins=bins, density=False)

    with open(output_file, 'w') as OUT:
        fieldnames = ['Bin start', 'Count']
        writer = csv.DictWriter(OUT, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for value, bin_start in zip(histogram[0], histogram[1]):
            writer.writerow({
                'Bin start': str(bin_start),
                'Count': str(value),
            })
