import math


#  TODO: move to utils
def read_chrom_sizes(input_file):

    chrom_sizes = dict()

    with open(input_file) as f:
        for line in f:
            chromosome, size = line.strip().split()
            chrom_sizes[chromosome] = int(size)

    return chrom_sizes


def write_corrected_bedgraph(input_bedgraph, chrom_sizes_fn, output_bedgraph,
                             y_int, scalar, mean_log, sd_log, slope):

    def print_line(chromosome, start, end, value, OUT):
        if value:  # Only prints non-zero values
            OUT.write(
                '{}\t{}\t{}\t{}\n'.format(
                    chromosome,
                    str(start - 1),
                    str(end),
                    str(value),
                )
            )

    chrom_sizes = read_chrom_sizes(chrom_sizes_fn)

    last_pos = None
    last_val = None
    last_chr = None
    start = None
    end = None

    with open(input_bedgraph) as f, open(output_bedgraph, 'w') as OUT:
        for line in f:
            chromosome, start, end, coverage = line.strip().split()[:4]
            start = int(start) + 1  # Convert from zero-based
            end = int(end)

            for position in range(start, end + 1):
                relative_pos = min(
                    position - 1, chrom_sizes[chromosome] - position)
                if relative_pos == 0:
                    relative_pos = 1

                corr = float(coverage) / (
                    scalar * (
                        math.erf((mean_log - math.log(relative_pos)) /
                        (math.sqrt(2) * sd_log))
                    ) + y_int + (slope * relative_pos)
                )
                value = int(round(corr, 0))

                if chromosome == last_chr and coverage == last_val and \
                        position == last_pos + 1:
                    end = position
                else:
                    print_line(last_chr, start, end, last_val, OUT)
                    start = position
                    end = position

                last_pos = position
                last_val = value
                last_chr = chromosome

        print_line(last_chr, start, end, last_val, OUT)
