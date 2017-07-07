import numpy
from scipy.stats import norm
from scipy.optimize import curve_fit

from fitting import gaussian


def calculate_depth_distribution(depths, output):

    depth_max = max(depths)
    hist = numpy.histogram(depths, bins=range(1, depth_max + 2), density=True)

    p0_mu, p0_sigma = norm.fit(depths)

    popt, pcov = curve_fit(gaussian, hist[1][:-1], hist[0], p0=[p0_mu, p0_sigma])
    mu, sigma = popt
    sigma = abs(sigma)

    # TODO: this is largely copied from calculate_bias_distribution
    with open(output, 'w') as OUT:

        OUT.write('Average depth: {}\n'.format(str(mu)))
        OUT.write('Standard deviation of depths: {}\n\n'.format(str(sigma)))

        OUT.write('Depth distribution:\n\n')
        OUT.write('\t'.join(['Depth', 'Frequency', 'Fit value']) +
            '\n')

        for hist_value, _bin in zip(hist[0], hist[1]):
            OUT.write('\t'.join((
                str(_bin),
                str(hist_value),
                str(norm.pdf(_bin, mu, sigma)),
            )) + '\n')

    return mu, sigma


def calculate_depth_distribution_bedgraph(in_bedgraph, output):

    depths = []

    with open(in_bedgraph) as f:
        for line in f:

            chromosome, start, end, coverage = line.strip().split()
            for i in range(int(start), int(end)):
                depths.append(int(coverage))

    return calculate_depth_distribution(depths, output)


def calculate_depth_distribution_mpileup(input_mpileup, output):

    depths = []

    with open(input_mpileup) as f:
        for line in f:

            line_split = line.strip().split()

            chromosome, position, reference_base, coverage = line_split[:4]
            position = int(position)
            coverage = int(coverage)
            if int(coverage) > 0:
                bases = line_split[4]
            else:
                bases = ''

            i = 0
            while i < len(bases):
                if bases[i] == '^':
                    i += 1
                elif bases[i] == '*':
                    coverage += 1
                i += 1

            if coverage > 0:
                depths.append(coverage)

    return calculate_depth_distribution(depths, output)


def process_chromosome_values(chromosome, chromosome_values, mu, sigma, OUT,
                              window=51, p_threshold=0.01):

    def write_position_to_filter(chromosome, position, depth, p):
        OUT.write('{}\t{}\t{}\n'.format(
            chromosome,
            str(i + 1),
            str(chromosome_values[i]),
            str(p),
        ))

    d = int((window - 1) / 2)
    norm_dist = norm(mu, sigma)

    keep_threshold = [mu, mu]
    filter_threshold = [float('-inf'), float('inf')]

    for i in range(d, len(chromosome_values) - d):

        window_depth = numpy.mean(chromosome_values[i - d:i + d + 1])

        if not (
            window_depth >= keep_threshold[0] and
            window_depth <= keep_threshold[1]
        ):
            if (
                window_depth <= filter_threshold[0] and
                window_depth >= filter_threshold[1]
            ):
                write_position_to_filter(
                    chromosome,
                    str(i + 1),
                    str(chromosome_values[i]),
                    p,
                )
            else:
                if window_depth < mu:
                    p = norm_dist.cdf(window_depth)

                    if p >= p_threshold:
                        keep_threshold[0] = window_depth
                    else:
                        filter_threshold[0] = window_depth
                        write_position_to_filter(
                            chromosome,
                            str(i + 1),
                            str(chromosome_values[i]),
                            p,
                        )

                elif window_depth > mu:
                    p = 1. - norm_dist.cdf(window_depth)

                    if p >= p_threshold:
                        keep_threshold[1] = window_depth
                    else:
                        filter_threshold[1] = window_depth
                        write_position_to_filter(
                            chromosome,
                            str(i + 1),
                            str(chromosome_values[i]),
                            p,
                        )


def filter_regions_by_depth(depths, chrom_sizes, mu, sigma,
                            filtered_regions_output):

    with open(filtered_regions_output, 'w') as OUT:

        for chromosome in sorted(depths.keys()):

            chromosome_depths = depths[chromosome]
            chromosome_values = []
            last_pos = 0

            for depth in chromosome_depths:

                for __ in range(last_pos + 1, depth['position']):
                    chromosome_values.append(0)

                chromosome_values.append(depth['coverage'])

                last_pos = depth['position']

            for __ in range(last_pos + 1, chrom_sizes[chromosome] + 1):
                chromosome_values.append(0)

            process_chromosome_values(
                chromosome, chromosome_values, mu, sigma, OUT)


def filter_regions_by_depth_bedgraph(bedgraph_file, chrom_sizes, mu,
                                     sigma, filtered_regions_output):
    depths = dict()

    with open(bedgraph_file) as f:
        for line in f:

            chromosome, start, end, coverage = line.strip().split()
            start = int(start) + 1  # Convert from zero-based
            end = int(end)
            coverage = int(coverage)

            if chromosome not in depths:
                depths[chromosome] = []

            for position in range(start, end + 1):
                depths[chromosome].append({
                    'position': position,
                    'coverage': coverage
                })

    filter_regions_by_depth(depths, chrom_sizes, mu, sigma,
        filtered_regions_output)


def filter_regions_by_depth_mpileup(mpileup_file, chrom_sizes, mu,
                                     sigma, filtered_regions_output):

    depths = dict()

    with open(mpileup_file) as f:
        for line in f:

            line_split = line.strip().split()

            chromosome, position, reference_base, coverage = line_split[:4]
            position = int(position)
            coverage = int(coverage)

            if chromosome not in depths:
                depths[chromosome] = []

            if int(coverage) > 0:
                bases = line_split[4]
            else:
                bases = ''

            i = 0
            while i < len(bases):
                if bases[i] == '^':
                    i += 1
                elif bases[i] == '*':
                    coverage += 1
                i += 1

            if coverage > 0:
                depths[chromosome].append({
                    'position': position,
                    'coverage': coverage,
                })

    filter_regions_by_depth(depths, chrom_sizes, mu, sigma,
        filtered_regions_output)
