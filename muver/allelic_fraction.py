#!/usr/bin/env python

import re
import math
import argparse
import subprocess

from collections import defaultdict

from wrappers.samtools import mpileup_iter


def get_allelic_fractions(bam_file, ref_fn, output_file):

    # #  Run samtools mpileup
    # proc = subprocess.Popen([
    #     samtools_path, 'mpileup',
    #     '-q', '5',
    #     '-Q', '10',
    #     '-B',
    #     '-d', '100000',
    #     '-f', ref_fn,
    #     bam_file,
    # ], stdout = subprocess.PIPE)
    #
    # fraction_histogram = defaultdict(int)
    #
    # for line in iter(proc.stdout.readline, ''):

    fraction_histogram = defaultdict(int)

    for line in mpileup_iter(bam_file, ref_fn):

        present_alleles = set()

        #  TODO: separate plus and minus tallies are unnecessary
        plus_tally = defaultdict(int)
        minus_tally = defaultdict(int)
        deletion_count = 0  # For '*' bases
        total = 0

        previous = None
        skip = None

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
                total += 1
            elif bases[i] == ',':
                present_alleles.add(reference)
                minus_tally[reference] += 1
                total += 1

            elif re.match('[ACGT]', bases[i]):
                present_alleles.add(bases[i])
                plus_tally[bases[i]] += 1
                total += 1
            elif re.match('[acgt]', bases[i]):
                present_alleles.add(bases[i].upper())
                minus_tally[bases[i]] += 1
                total += 1

            elif re.match('[+-]', bases[i]):
                indel_type = bases[i]
                i += 1

                indel_length = int(bases[i])
                i += 1

                indel = indel_type + bases[i:i+indel_length].upper()
                present_alleles.add(indel)
                if re.match('[ACGT]', bases[i:i+indel_length]):
                    plus_tally[indel] += 1
                elif re.match('[acgt]', bases[i:i+indel_length]):
                    minus_tally[indel] += 1
                total += 1
                i += indel_length

            elif bases[i] == '^':
                i += 1
            elif bases[i] == '*':
                deletion_count += 1
                total += 1

            i += 1

        if total > 20:
            for allele in present_alleles:
                fraction = float(plus_tally[allele] + minus_tally[allele]) / total
                fraction_histogram[str(round(fraction, 2))] += 1
            if deletion_count > 0:
                fraction = float(deletion_count) / total
                fraction_histogram[str(round(fraction, 2))] += 1

    with open(output_file, 'w') as OUT:
        for i in range(101):
            display = str(round(float(i)/100, 2))
            OUT.write('{}\t{}\n'.format(
                display,
                fraction_histogram[display],
            ))
