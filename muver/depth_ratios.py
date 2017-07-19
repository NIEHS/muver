import csv
import math
import numpy
import os
from collections import defaultdict

from depth_distribution import calculate_depth_distribution_bedgraph
from reference import read_chrom_sizes


def calculate_depth_ratios(input_bedgraph, reference_assembly, output_file,
                           mean=None, bin_size=500):
    '''
    For binned depths across the input bedGraph file, calculate the ratio
    relative to the mean.
    '''
    if not mean:
        mean = calculate_depth_distribution_bedgraph(
            input_bedgraph, os.devnull)[0]
    chrom_sizes = read_chrom_sizes(reference_assembly)

    binned_values = defaultdict(list)

    with open(input_bedgraph) as f:
        for line in f:
            chromosome, start, end, depth = line.strip().split('\t')
            start = int(start)
            end = int(end)
            depth = float(depth)

            if depth != 0:
                for position in range(start + 1, end + 1):
                    distance_from_chrom_start = end - 1
                    distance_from_chrom_end = chrom_sizes[chromosome] - end
                    min_distance = min(distance_from_chrom_start,
                                   distance_from_chrom_end)

                    if depth / mean >= 0.25 and depth / mean <= 4:
                        bin_index = math.floor(float(min_distance) / 500)
                        _bin = (
                            bin_index * 500,
                            (bin_index + 1) * 500 - 1,
                        )
                        binned_values[_bin].append(depth / mean)

    output_values = []
    for _bin, values in sorted(binned_values.items(), key=lambda x: x[0][0]):
        output_values.append({
            'Bin Start': int(_bin[0]),
            'Bin End': int(_bin[1]),
            'Bin Center': int(math.ceil(float(sum(_bin)) / 2)),
            'Median Value': numpy.median(values),
        })

    with open(output_file, 'w') as OUT:
        fieldnames = ['Bin Start', 'Bin End', 'Bin Center', 'Median Value']
        writer = csv.DictWriter(OUT, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in output_values:
            writer.writerow(row)
