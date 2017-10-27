from collections import defaultdict
import csv
import matplotlib.pyplot as plt
import math
import re
from scipy.optimize import curve_fit

from fitting import logistic
from wrappers.samtools import view_bam


def calculate_repeat_occurrences(repeats):
    '''
    Iterate through a repeats dictionary. Count the number of occurrences of
    repeats of a given repeat unit length and repeat tract length. Return as
    dict.

    Only count full occurrences of a given repeat unit, ignoring trailing
    bases.

    e.g. GAGAGAGAG -- repeat_unit_length: 2, repeat_tract_length: 8
    '''
    occurrences = defaultdict(int)

    for chromosome_repeats in repeats.values():
        for position, position_repeats in chromosome_repeats.items():
            for repeat in position_repeats:

                if position == repeat['start']:

                    repeat_unit_length = len(repeat['unit'])
                    repeat_tract_length = repeat_unit_length * \
                        repeat['sequence'].count(repeat['unit'])

                    occurrences[(repeat_unit_length, repeat_tract_length)] += 1

    return occurrences


def calculate_repeat_indel_counts(repeats, sam_iter):
    '''
    Iterating through a BAM/SAM file, count the number of indels with repeats,
    noting the repeat unit length and repeat tract length. Return as dict
    with overall depth counts and insertion and deletion counts.
    '''
    counts = dict()
    for field in ['depth', 'insertion', 'deletion']:
        counts[field] = dict()
        for repeat_length in range(1, 5):
            counts[field][repeat_length] = defaultdict(int)

    for line in sam_iter:

        if line.startswith('@'):
            pass

        else:

            line_split = line.strip().split('\t')

            chromosome = line_split[2]
            position = int(line_split[3])
            cigar_string = line_split[5]
            sequence = line_split[9]

            total_length = 0
            length = 0
            _sum = 0

            starts = []
            deletions = []
            lengths = []
            ops = []

            for match in re.finditer('(\d+)([MISD])', cigar_string):

                try:
                    starts.append(starts[-1] + length)
                except IndexError:
                    starts.append(length)

                length = int(match.group(1))
                op = match.group(2)

                ops.append(op)

                if op == 'M' or op == 'D':
                    total_length += length

                if op == 'D':
                    deletions.append(length)
                    length = 0
                else:
                    deletions.append(0)

                lengths.append(length)

            end = position + total_length - 1

            for i, op in enumerate(ops):

                if op == 'S':
                    _sum += lengths[i]

                elif op == 'I':
                    left = position + starts[i] - 1
                    left -= _sum

                    try:
                        position_repeats = repeats[chromosome][left]
                    except KeyError:
                        position_repeats = []

                    for repeat in position_repeats:
                        if repeat['start'] == left:

                            if repeat['start'] + len(repeat['sequence']) < end:

                                repeat_unit_length = len(repeat['unit'])
                                repeat_length = repeat_unit_length * \
                                    repeat['sequence'].count(repeat['unit'])

                                counts['depth'][repeat_unit_length][repeat_length] += 1
                                counts['insertion'][repeat_unit_length][repeat_length] += 1

                    _sum += lengths[i]

                elif op == 'D':
                    left = position + starts[i] - 1
                    left -= _sum

                    try:
                        position_repeats = repeats[chromosome][left]
                    except KeyError:
                        position_repeats = []

                    for repeat in position_repeats:
                        if repeat['start'] == left:

                            if repeat['start'] + len(repeat['sequence']) < end:

                                repeat_unit_length = len(repeat['unit'])
                                repeat_length = repeat_unit_length * \
                                    repeat['sequence'].count(repeat['unit'])

                                counts['depth'][repeat_unit_length][repeat_length] += 1
                                counts['deletion'][repeat_unit_length][repeat_length] += 1

                    _sum -= deletions[i]

                elif op == 'M':
                    for j in range(lengths[i]):
                        left = position + starts[i] - _sum + j

                        try:
                            position_repeats = repeats[chromosome][left]
                        except KeyError:
                            position_repeats = []

                        for repeat in position_repeats:
                            if repeat['start'] == left and j < lengths[i] - 1:

                                if repeat['start'] + len(repeat['sequence']) < end:

                                    repeat_unit_length = len(repeat['unit'])
                                    repeat_length = repeat_unit_length * \
                                        repeat['sequence'].count(repeat['unit'])

                                    counts['depth'][repeat_unit_length][repeat_length] += 1

    return counts


def calculate_repeat_indel_rates(indel_counts, repeat_occurrences,
                                 occurrence_filter=10):
    '''
    Considering overall depth counts, find the rates for insertions and
    deletions within repeats. Store values as a function of repeat length and
    repeat tract length.

    Do not report a rate if repeats of a given repeat unit length and repeat
    tract length appear fewer times than the occurrence filter. Read these
    occurrences from the repeat_occurrences dict.
    '''
    rates = dict()
    for field in ['insertion', 'deletion']:
        rates[field] = dict()
        for repeat_length in range(1, 5):
            rates[field][repeat_length] = dict()

    for event in ['insertion', 'deletion']:
        for repeat_length, counts in indel_counts[event].items():
            for tract_length, count in counts.items():

                if repeat_occurrences[(repeat_length, tract_length)] >= \
                        occurrence_filter:

                    a = indel_counts[event][repeat_length][tract_length]
                    b = indel_counts['depth'][repeat_length][tract_length]

                    rates[event][repeat_length][tract_length] = float(a) / b

    return rates


def fit_rates(indel_rates):
    '''
    Considering repeats of a given repeat unit length, fit rates to a logistic
    function of repeat tract length.  Return fit parameters in a dict.
    '''
    fits = dict()

    for event in ['insertion', 'deletion']:
        fits[event] = dict()

        for repeat_length, rates in indel_rates[event].items():

            tract_lengths = []
            repeat_rates = []

            for tract_length, rate in rates.items():

                tract_lengths.append(tract_length)
                repeat_rates.append(math.log10(rate))

            if repeat_rates:

                p0_k = 0.1
                p0_L = max(repeat_rates) - min(repeat_rates)
                p0_M = min(repeat_rates)

                mid_value = (max(repeat_rates) + min(repeat_rates)) / 2
                mid_diff = float('inf')
                for i, val in enumerate(repeat_rates):
                    diff = abs(val - mid_value)
                    if diff < mid_diff:
                        mid_diff = diff
                        p0_x0 = tract_lengths[i]

                popt, pcov = curve_fit(
                    logistic,
                    tract_lengths,
                    repeat_rates,
                    p0=[p0_x0, p0_L, p0_M, p0_k],
                    max_nfev=1000000,
                    method='trf',
                )
                x0, L, M, k = popt

                fits[event][repeat_length] = {
                    'x0': x0,
                    'L': L,
                    'M': M,
                    'k': k,
                }

    return fits


def plot_fits(indel_rates, fits, output_header):
    '''
    Considering logistic function fits, compare fit values to observed rates.
    Fit values and observed rates are reported in scatter plots written to PNG
    files named using output_header.
    '''
    for event, event_fits in fits.items():
        for repeat_length, repeat_fit in event_fits.items():

            rates = indel_rates[event][repeat_length]

            tract_lengths = []
            fitted_values = []
            raw_values = []

            for x in sorted(rates.keys()):
                tract_lengths.append(x)
                raw_values.append(math.log10(rates[x]))

                fitted = logistic(
                    x,
                    repeat_fit['x0'],
                    repeat_fit['L'],
                    repeat_fit['M'],
                    repeat_fit['k'],
                )
                fitted_values.append(fitted)

            plt.plot(sorted(tract_lengths), raw_values, 'ko',
                     sorted(tract_lengths), fitted_values, 'b')
            plt.savefig('{}_{}_{}.png'.format(
                output_header, event, str(repeat_length)))
            plt.close()


def print_fits(fits, output_file):
    '''
    Print parameters from logistic function fits in the input dict to an output
    file.
    '''
    ordered_fits = []

    for event in ('insertion', 'deletion'):
        for repeat_length in (1, 2, 3, 4):

            if repeat_length in fits[event]:

                ordered_fits.append({
                    'Event': event,
                    'Repeat length': str(repeat_length),
                })
                ordered_fits[-1].update(fits[event][repeat_length])

    field_names = ('Event', 'Repeat length', 'x0', 'L', 'M', 'k')

    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=field_names, delimiter='\t')

        writer.writeheader()
        for fit in ordered_fits:
            writer.writerow(fit)

def print_rates(indel_rates, output_file):
    '''
    Print observed indel error rates to an output file.
    '''
    with open(output_file.strip('.txt') + '.rates.txt','w') as f:

        max = 0
        for repeat_length in (1, 2, 3, 4):
            for event in ('insertion', 'deletion'):
                f.write('\t' + event + '_' + str(repeat_length))
                curmax = sorted(indel_rates[event][repeat_length].keys(), \
                    reverse=True)[0]
                if curmax > max:
                    max = curmax

        f.write('\n')

        i = 4
        while i <= max:
            count = 0
            out = str(i)
            for repeat_length in (1, 2, 3, 4):
                for event in ('insertion', 'deletion'):
                    if i in indel_rates[event][repeat_length]:
                        rate = indel_rates[event][repeat_length][i]
                        out += '\t'
                        out += str(math.log10(rate))
                        count = count + 1
                    else:
                        out += '\tX'
            if count > 0:
                f.write(out + '\n')
            i = i + 1

def read_fits(fits_file):
    '''
    Read logistic function parameters from a tab-delimited TXT file.
    '''
    fits = dict()

    with open(fits_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:

            if row['Event'] not in fits:
                fits[row['Event']] = dict()

            fits[row['Event']][int(row['Repeat length'])] = {
                'x0': float(row['x0']),
                'L': float(row['L']),
                'M': float(row['M']),
                'k': float(row['k']),
            }

    return fits

def read_rates(rates_file):
    '''
    Read observed indel error rates from a tab-delimited TXT file.
    '''
    rates=dict()

    for t in ['insertion', 'deletion']:
        rates[t]=dict()
        for n in range(1,5):
            rates[t][n]=dict()

    with open(rates_file) as f:
        f.readline()
        for line in f:
            fields=line.strip().split('\t')
            i=1
            for n in range(1,5):
                for t in ['insertion', 'deletion']:
                    if fields[i] != 'X':
                        rates[t][n][int(fields[0])] = 10**float(fields[i])
                    i += 1
    return(rates)

def fit_repeat_indel_rates(repeats, bam_file, output_file,
                                    output_plot_header=None):
    '''
    Considering a dictionary of repeats, iterate through a BAM/SAM file and
    record occurrences of indels within repeat sequences. Then, find rates
    of indel occurrences as a function of repeat unit length and repeat tract
    length. For each repeat unit length, fit observed rates as a function of
    repeat tract length to a logistic function. Write the parameters of these
    fits to an output file. If output_plot_header is specified, plot fitted
    values for visual validation.
    '''
    bam_iter = view_bam(bam_file)

    repeat_occurrences = calculate_repeat_occurrences(repeats)

    indel_counts = calculate_repeat_indel_counts(repeats, bam_iter)
    indel_rates = calculate_repeat_indel_rates(
        indel_counts, repeat_occurrences)

    fits = fit_rates(indel_rates)

    print_fits(fits, output_file)
    print_rates(indel_rates, output_file)
    if output_plot_header:
        plot_fits(indel_rates, fits, output_plot_header)
