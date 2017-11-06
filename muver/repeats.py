from collections import OrderedDict
import copy
import itertools
import regex as re
import random


def reverse_enumerate(iterable):
    '''
    Enumerate through an iterable in reverse, reporting the index consistent
    with the original iterable.
    '''
    return itertools.izip(reversed(xrange(len(iterable))), reversed(iterable))


def generate_repeat_units():
    '''
    Given canonical bases, generate a set of all possible repeat units up to
    a length of four.
    '''
    bases = ['A', 'C', 'G', 'T']
    repeat_units = set(copy.copy(bases))

    for i in range(3):
        for unit in copy.copy(repeat_units):
            for base in bases:
                repeat_units.add(unit + base)

    temp_repeat_units = set()
    for repeat_unit in repeat_units:
        keep = True

        n = len(repeat_unit)
        if n > 1:
            if repeat_unit[0] * (n - 1) == repeat_unit[1:n]:
                keep = False
        if n == 4:
            if repeat_unit[0:2] == repeat_unit[2:4]:
                keep = False

        if keep:
            temp_repeat_units.add(repeat_unit)
    repeat_units = temp_repeat_units

    return repeat_units


def check_repeats(repeat_1, repeat_2):
    '''
    Check to see if repeat_1 is a possible permutation of repeat_2.

    e.g. check_repeats('AGCT', 'GCTA') is True, check_repeats('AGCT', 'ATGC')
    is False.
    '''
    if repeat_1 == repeat_2:

        return True

    elif len(repeat_1) == len(repeat_2):

        for i in range(1, len(repeat_1)):

            shuffled_repeat = repeat_1[i:] + repeat_1[:i]

            if shuffled_repeat == repeat_2:
                return True

    return False


def create_repeat_file(fasta_file, output_file):
    '''
    For a given FASTA file, enumerate all repeats to an output file.
    '''
    repeat_units = generate_repeat_units()
    sequences = OrderedDict()
    seq_name = None
    seq = ''

    groups = dict()

    for repeat_unit in repeat_units:

        groups[repeat_unit] = dict()

        for other_repeat_unit in repeat_units:

            groups[repeat_unit][other_repeat_unit] = \
                check_repeats(repeat_unit, other_repeat_unit)

    with open(fasta_file) as f:

        for line in f:

            if line.startswith('>'):  # New FASTA entry

                sequences[seq_name] = seq
                seq = ''
                seq_name = ''.join(line.split('>')[1:]).strip()
                sequences[seq_name] = ''

            else:

                seq += line.strip()

    sequences[seq_name] = seq

    with open(output_file, 'w') as OUT:

        for sequence_name, sequence in sequences.items():

            matches = []

            for repeat_unit in repeat_units:

                repeat_length = len(repeat_unit)

                fragments = []
                for i in range(1, len(repeat_unit)):
                    fragments.append(repeat_unit[:-i])

                search_pattern = '({}){{2,}}({}){{0,1}}'.format(
                    repeat_unit,
                    '|'.join(fragments),
                )

                last_start = None
                for match in re.finditer(search_pattern, sequence,
                        overlapped=True):

                    keep = True

                    if last_start:
                        if match.start() - repeat_length == last_start:
                            keep = False

                    if keep:
                        matches.append({
                            'sequence': match.group(0),
                            'repeat_unit': repeat_unit,
                            'start': match.start(),
                            'end': match.end(),
                        })

                    last_start = match.start()

            sort = sorted(matches, key=lambda x: (x['start'], -x['end']))
            kept_matches = []

            i = len(sort) - 1
            while i >= 0:

                keep = True

                j = i - 1
                while j >= 0:

                    if (
                        sort[i]['start'] >= sort[j]['start'] and
                        sort[i]['end'] <= sort[j]['end'] and
                        groups[sort[i]['repeat_unit']][sort[j]['repeat_unit']]
                    ):
                        keep = False
                        break

                    if sort[i]['start'] > sort[j]['end']:
                        break

                    j = j - 1

                if keep:
                    kept_matches.append(sort[i])

                i = i - 1

            for match in sorted(kept_matches, key=lambda x: x['start']):

                OUT.write('\t'.join((
                    sequence_name,
                    match['sequence'],
                    str(len(match['repeat_unit'])),
                    match['repeat_unit'],
                    str(match['start']),
                    str(match['end']),
                )) + '\n')


def extract_repeat_file_sample(repeat_file, sample_file, total):
    '''
    Extract a random sample of repeat loci from a genome-wide list
    '''
    with open(repeat_file, 'r', 1) as f:
        for i, l in enumerate(f):
            pass
        i += 1

    keep = dict(zip(random.sample(range(0, i), total),itertools.repeat(0)))

    with open(repeat_file, 'r', 1) as f, open(sample_file, 'w') as OUT:
        for x, line in enumerate(f):
            if x in keep:
                OUT.write(line)
