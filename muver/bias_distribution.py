from collections import defaultdict
import math
import numpy
import re
from scipy.stats import norm
from scipy.optimize import curve_fit

from fitting import gaussian
from wrappers import samtools


def calculate_bias_distribution(_iter, ref_fn, output):
    '''
    For an iterable, create a histogram of strand bias values.  Then
    fit those values using a log-normal distribution. Output parameters
    for the fit and the histogram in a TXT file.

    _iter may be generated using calculate_bias_distribution_bam or
    calculate_bias_distribution_mpileup.
    '''
    log_ratios = []

    for line in _iter:
        present_alleles = set()

        plus_tally = defaultdict(int)
        minus_tally = defaultdict(int)

        line_split = line.strip().split()

        reference, coverage = line_split[2:4]
        if int(coverage) > 0:
            bases = line_split[4]
        else:
            bases = ''

        i = 0
        while i < len(bases):

            if bases[i] == '.':
                present_alleles.add(reference)
                plus_tally[reference] += 1
            elif bases[i] == ',':
                present_alleles.add(reference)
                minus_tally[reference] += 1

            elif re.match('[ACGTN]', bases[i]):
                present_alleles.add(bases[i])
                plus_tally[bases[i]] += 1
            elif re.match('[acgtn]', bases[i]):
                present_alleles.add(bases[i].upper())
                minus_tally[bases[i].upper()] += 1

            elif re.match('[+-]', bases[i]):
                indel_type = bases[i]
                i += 1

                num = ''
                while re.match('[0-9]', bases[i]):
                    num += bases[i]
                    i += 1

                indel_length = int(num)

                indel = indel_type + bases[i:i+indel_length].upper()
                present_alleles.add(indel)
                if re.match('[ACGTN]', bases[i:i+indel_length]):
                    plus_tally[indel] += 1
                elif re.match('[acgtn]', bases[i:i+indel_length]):
                    minus_tally[indel] += 1
                i += indel_length - 1

            elif bases[i] == '^':
                i += 1

            i += 1

        for allele in present_alleles:
            if allele in plus_tally and allele in minus_tally:
                ratio = float(plus_tally[allele]) / minus_tally[allele]
                log_ratios.append(math.log(ratio))

    p0_mu, p0_sigma = norm.fit(log_ratios)
    if p0_sigma == 0:
        p0_sigma = 0.01

    hist = numpy.histogram(log_ratios,
        bins=[float(x)/10 for x in range(-50, 51)], density=True)
    popt, pcov = curve_fit(gaussian, hist[1][:-1], hist[0],
        p0=[p0_mu, p0_sigma], maxfev=100000)
    mu, sigma = popt
    sigma = abs(sigma)

    with open(output, 'w') as OUT:

        OUT.write('Average log ratio: {}\n'.format(str(mu)))
        OUT.write(
            'Standard deviation of log ratios: {}\n\n'.format(str(sigma)))

        OUT.write('Bias distribution:\n\n')
        OUT.write('\t'.join(['Strand log ratio', 'Frequency', 'Fit value']) +
            '\n')

        for hist_value, _bin in zip(hist[0], hist[1]):
            OUT.write('\t'.join((
                str(round(_bin, 1)),
                str(hist_value),
                str(norm.pdf(_bin, mu, sigma)),
            )) + '\n')

    return mu, sigma


def calculate_bias_distribution_bam(input_bam, ref_fn, output):
    '''
    Perform bias distribution characterization with a BAM file.
    '''
    return calculate_bias_distribution(
        samtools.mpileup_iter(input_bam, ref_fn), ref_fn, output)


def calculate_bias_distribution_mpileup(input_mpileup, ref_fn, output):
    '''
    Perform bias distribution characterization with an mpileup TXT file.
    '''
    with open(input_mpileup) as f:

        return calculate_bias_distribution(
            iter(f.readline, ''), ref_fn, output)
