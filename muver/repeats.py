import copy
import argparse
import regex as re
import itertools

from collections import OrderedDict


def reverse_enumerate(iterable):
    return itertools.izip(reversed(xrange(len(iterable))), reversed(iterable))


def generate_repeat_units():

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

    if repeat_1 == repeat_2:

        return True

    elif len(repeat_1) == len(repeat_2):

        for i in range(1, len(repeat_1)):

            shuffled_repeat = repeat_1[i:] + repeat_1[:i]

            if shuffled_repeat == repeat_2:
                return True

    return False


def create_repeat_file(fasta_file, output_file):

    repeat_units = generate_repeat_units()
    sequences = OrderedDict()
    seq_name = None

    with open(fasta_file) as f:

        for line in f:

            if line.startswith('>'):  # New FASTA entry

                seq_name = ''.join(line.split('>')[1:]).strip()
                sequences[seq_name] = ''

            else:

                sequences[seq_name] += line.strip()

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

            for i, match in reverse_enumerate(sort):

                keep = True

                for other_match in reversed(sort[:i]):

                    if match['start'] >= other_match['start'] and \
                            match['end'] <= other_match['end'] and \
                            check_repeats(match['repeat_unit'], other_match['repeat_unit']):

                        keep = False
                        break

                    if match['start'] > other_match['end']:
                        break

                if keep:
                    kept_matches.append(match)

            for match in sorted(kept_matches, key=lambda x: x['start']):

                OUT.write('\t'.join((
                    sequence_name,
                    match['sequence'],
                    str(len(match['repeat_unit'])),
                    match['repeat_unit'],
                    str(match['start']),
                    str(match['end']),
                )) + '\n')
