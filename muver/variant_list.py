import math
import numpy

from variant import Variant
from utils import read_excluded_regions, read_repeats_var, read_filtered_sites


def get_allele_values(alleles, in_dict):
    '''
    For a nested dict, where the first level has keys by allele and the second
    level has keys by strand, return a flattened list of those values
    respecting the order of the input allele list.
    '''
    d = []
    for allele in alleles:
        for strand in ['forward', 'reverse']:
            d.append(in_dict[allele][strand])
    return d


def format_value(obj):
    '''
    Coerce the input object into a format amenable to legible printing.
    '''
    if obj is None:
        return 'NA'
    elif type(obj) == bool:
        return str(int(obj))
    elif type(obj) == list or type(obj) == tuple:
        return ','.join([str(x) for x in obj])
    elif type(obj) == int or type(obj) == float or type(obj) == numpy.float64:
        return str(obj)
    elif type(obj) == str:
        return obj


def get_vcf_genotype(genotype, alleles, ploidy):
    '''
    Given the genotype, allele list, and ploidy, return the genotype in VCF
    format.
    '''
    if genotype:
        allele_indices = [str(alleles.index(a)) for a in genotype]
        return '/'.join(allele_indices)
    else:
        return '/'.join(['.'] * ploidy)


def get_vcf_mutations(mutations):
    '''
    Format a list of mutations for writing to a VCF file.
    '''
    if mutations:
        names = [m['name'] for m in mutations]
        return ','.join(names)
    else:
        return '.'


def format_vcf_field(field):
    '''
    Coerce objects to a format for writing to a VCF file.
    '''
    if field is None:
        return '.'
    elif type(field) == bool:
        return str(int(field))
    elif type(field) == list or type(field) == tuple:
        return ','.join([str(x) for x in field])
    else:
        return str(field)


def read_vcf_filter(flags):
    '''
    Given input flags, return a VCF filter field.
    '''
    filters = []
    for flag, flag_name in flags:
        if flag:
            filters.append(flag_name)

    if filters:
        return ';'.join(filters)
    else:
        return 'PASS'


def write_vcf_info(f_field, f_id, f_number, f_type, f_description, OUT):
    '''
    Write a VCF header line to the OUT file. f_field, f_id, f_number, and
    f_type correspond to fields specified in VCF documentation.
    '''
    fields = []
    for field, value in (
        ('ID', f_id),
        ('Number', f_number),
        ('Type', f_type),
        ('Description', f_description),
    ):
        if value is not None:
            if field == 'Description':
                value = '\"{}\"'.format(value)
            fields.append('{}={}'.format(field, value))

    OUT.write('##{}=<{}>\n'.format(f_field, ','.join(fields)))


class VariantList(object):
    '''
    In addition to containing a list of variants, also contains other
    information for a MuVer run.
    '''
    def __init__(self, input_vcf_fn, samples, excluded_regions_fn, repeats_fn,
                 control_sample, chrom_sizes, depth_threshold=20):

        self.variants = []
        self.input_vcf_fn = input_vcf_fn

        self.samples = samples
        self.control_sample = control_sample
        self.excluded_regions_fn = excluded_regions_fn
        self.repeats_fn = repeats_fn

        genome_size = sum(chrom_sizes.values())
        self.p_threshold = 1.0 - ((1.0 - 0.01) ** (1.0 / genome_size))

        self.depth_threshold = depth_threshold

        if self.excluded_regions_fn:
            excluded_regions = read_excluded_regions(self.excluded_regions_fn)
        else:
            excluded_regions = set()

        filtered_sites = read_filtered_sites(self.samples)

        self.read_variants_from_vcf()

        var_dict=dict()

        for variant in self.variants:
            var_dict[(variant.chromosome, variant.position)] = 0

        repeats = read_repeats_var(self.repeats_fn, var_dict)

        for variant in self.variants:

            variant.check_if_excluded(excluded_regions)
            variant.check_depth_threshold(self.depth_threshold)
            variant.check_filtered_sites(filtered_sites)
            variant.check_allele_coverage()
            variant.binomial_test()
            variant.chisquare_test()
            variant.check_if_significant(self.p_threshold)
            variant.assign_ploidy()
            variant.intersect_with_repeats(repeats)
            variant.find_repeat_expanded_alleles()
            variant.call_genotypes_and_subclonal_alleles()
            variant.subclonal_strand_bias_log_normal_test(samples)
            variant.subclonal_binomial_test(self.p_threshold)
            variant.subclonal_strand_bias_binomial_test()
            variant.call_mutations()
            variant.set_report_mutations_flag(self.p_threshold)
            variant.check_subclonal_validity(self.p_threshold)

    def __iter__(self):
        return iter(self.variants)

    def read_variants_from_vcf(self):
        '''
        Iterates through a GATK HaplotypeCaller VCF and creates a variant for
        each entry.
        '''
        with open(self.input_vcf_fn) as f:
            for line in f:

                if line.startswith('##'):
                    pass

                elif line.startswith('#'):
                    sample_names = line.strip().split('\t')[9:]
                    if self.control_sample.sample_name not in sample_names:
                        raise ValueError('Control sample not in input VCF.')

                    sample_indices = dict()
                    for i, sample_name in enumerate(sample_names):
                        sample_indices[i] = next(
                            (x for x in self.samples
                                if sample_name == x.sample_name),
                            None,
                        )

                else:
                    (chromosome, position, _id, ref, alt, qual, _filter, info,
                        _format) = line.strip().split('\t')[:9]

                    sample_fields = line.strip().split('\t')[9:]
                    alleles = [ref] + alt.split(',')
                    _format = _format.split(':')

                    sac = dict()
                    for sample in self.samples:
                        sac[sample] = dict()

                    try:
                        sac_index = _format.index('SAC')
                    except ValueError:
                        sac_index = None

                    for i, fields in enumerate(sample_fields):
                        sample = sample_indices[i]
                        if sac_index == None:
                            counts = [0] * len(alleles) * 2
                        else:
                            try:
                                counts = fields.split(':')[sac_index].split(',')
                            except IndexError:
                                counts = [0] * len(alleles) * 2

                        for allele in alleles:
                            forward = int(counts.pop(0))
                            reverse = int(counts.pop(0))
                            sac[sample][allele] = {
                                'forward': forward,
                                'reverse': reverse,
                            }

                    self.variants.append(Variant(chromosome, position, alleles,
                        self.samples, ref, self.control_sample, sac))

    def write_output_table(self, output_fn):
        '''
        Writes a line to an output file for each variant in the variant list.
        '''
        output_headers = [
            'Chromosome',
            'Position',
            'Excluded Region',
            'Alleles',
            'Intersects Repeat',
            'Repeat Correction Applied',
            '{} Depths'.format(self.control_sample.sample_name),
            '{} Allele Count Flag'.format(self.control_sample.sample_name),
            '{} Depth Flag'.format(self.control_sample.sample_name),
            '{} Filter Flag'.format(self.control_sample.sample_name),
            '{} Ploidy'.format(self.control_sample.sample_name),
            '{} Genotype'.format(self.control_sample.sample_name),
            '{} Subclonal Genotype'.format(self.control_sample.sample_name),
            '{} Subclonal Frequency'.format(self.control_sample.sample_name),
            '{} Genotyping Score'.format(self.control_sample.sample_name),
            '{} Subclonal Valid Flag'.format(self.control_sample.sample_name),
        ]

        for sample in [s for s in self.samples if s != self.control_sample]:
            output_headers.extend([
                '{} Depths'.format(sample.sample_name),
                '{} Allele Count Flag'.format(sample.sample_name),
                '{} Depth Flag'.format(sample.sample_name),
                '{} Filter Flag'.format(sample.sample_name),
                '{} Composite Score'.format(sample.sample_name),
                '{} Read Difference Flag'.format(sample.sample_name),
                '{} Ploidy'.format(sample.sample_name),
                '{} Genotype'.format(sample.sample_name),
                '{} Subclonal Genotype'.format(sample.sample_name),
                '{} Subclonal Frequency'.format(sample.sample_name),
                '{} Genotyping Score'.format(sample.sample_name),
                '{} Subclonal Valid Flag'.format(sample.sample_name),
                '{} Mutations'.format(sample.sample_name),
                '{} PAC Flag'.format(sample.sample_name),
            ])

        with open(output_fn, 'w') as OUTPUT:

            OUTPUT.write('\t'.join(output_headers) + '\n')

            for v in self.variants:

                sac = v.strand_allele_counts

                variant_fields = [
                    v.chromosome,
                    v.position,
                    v.excluded_flag,
                    v.alleles,
                    bool(v.intersected_repeats),
                    bool(v.intersected_repeat_added_sequence),
                    get_allele_values(v.alleles, sac[self.control_sample]),
                    v.allele_coverage_flags[self.control_sample],
                    v.depth_flags[self.control_sample],
                    v.filtered_sites_flags[self.control_sample],
                    v.sample_ploidy[self.control_sample],
                    v.sample_genotypes[self.control_sample],
                    v.sample_subclonals[self.control_sample]['genotype'],
                    v.sample_subclonals[self.control_sample]['frequency'],
                    v.sample_genotype_min_log_ratio_sum[self.control_sample],
                    v.sample_subclonal_valid_flag[self.control_sample],
                ]

                for sample in [s for s in self.samples if s != self.control_sample]:

                    muts = v.sample_called_mutations[sample]
                    if muts and v.report_mutations[sample]:
                        sample_mutations = [m['name'] for m in muts]
                        sample_called_loh = v.sample_called_loh[sample]
                    else:
                        sample_mutations = None
                        sample_called_loh = None

                    sample_fields = [
                        get_allele_values(v.alleles, sac[sample]),
                        v.allele_coverage_flags[sample],
                        v.depth_flags[sample],
                        v.filtered_sites_flags[sample],
                        v.sample_composite_p_values[sample],
                        v.sample_significance_flags[sample],
                        v.sample_ploidy[sample],
                        v.sample_genotypes[sample],
                        v.sample_subclonals[sample]['genotype'],
                        v.sample_subclonals[sample]['frequency'],
                        v.sample_genotype_min_log_ratio_sum[sample],
                        v.sample_subclonal_valid_flag[sample],
                        sample_mutations,
                        sample_called_loh,
                    ]

                    variant_fields.extend(sample_fields)

                OUTPUT.write('\t'.join([format_value(x) for x in variant_fields]) + '\n')

    def write_output_vcf(self, output_fn):
        '''
        Writes a VCF file describing variants in the list.
        '''
        samples = [self.control_sample] + \
            [s for s in self.samples if s != self.control_sample]

        with open(output_fn, 'w') as OUTPUT:

            write_vcf_info('INFO', 'IntersectsRepeat', '0', 'Flag', 'Intersects with repeat', OUTPUT)
            write_vcf_info('INFO', 'RepeatCorrectionApplied', '0', 'Flag', 'Repeat correction applied', OUTPUT)

            write_vcf_info('FILTER', 'ExcRegion', None, None, 'Intersects with excluded region', OUTPUT)

            write_vcf_info('FORMAT', 'SAC', '.', 'Integer', 'Strand allele counts', OUTPUT)
            write_vcf_info('FORMAT', 'ACF', '0', 'Flag', 'Allele count flag', OUTPUT)
            write_vcf_info('FORMAT', 'DF', '0', 'Flag', 'Depth flag', OUTPUT)
            write_vcf_info('FORMAT', 'PD', '1', 'Integer', 'Ploidy, specified by user', OUTPUT)
            write_vcf_info('FORMAT', 'GT', '1', 'String', 'Genotype', OUTPUT)
            write_vcf_info('FORMAT', 'SGT', '1', 'String', 'Subclonal genotype', OUTPUT)
            write_vcf_info('FORMAT', 'SF', '1', 'Float', 'Subclonal frequency', OUTPUT)
            write_vcf_info('FORMAT', 'SV', '0', 'Flag', 'Subclonal valid flag', OUTPUT)
            write_vcf_info('FORMAT', 'RD', '0', 'Flag', 'Read difference flag', OUTPUT)
            write_vcf_info('FORMAT', 'MT', '.', 'String', 'Called mutations', OUTPUT)
            write_vcf_info('FORMAT', 'PAC', '.', 'String', 'Called PAC flags', OUTPUT)

            headers = [
                'CHROM',
                'POS',
                'ID',
                'REF',
                'ALT',
                'QUAL',
                'FILTER',
                'INFO',
                'FORMAT',
            ]
            headers += [s.sample_name for s in samples]
            OUTPUT.write('#' + '\t'.join(headers) + '\n')

            for variant in self.variants:

                info_fields = []

                filter_fields = read_vcf_filter([
                    (variant.excluded_flag, 'ExcRegion'),
                ])

                info_fields.append('='.join((
                    'IntersectsRepeat', format_vcf_field(variant.intersect_repeat_flag))))
                info_fields.append('='.join((
                    'RepeatCorrectionApplied', format_vcf_field(bool(variant.intersected_repeat_unit)))))

                info_fields = ';'.join(info_fields)

                alt_alleles = \
                    [a for a in variant.alleles if a != variant.ref_allele]
                alleles = [variant.ref_allele] + alt_alleles

                format_fields = [
                    'SAC', 'ACF', 'DF', 'PD', 'GT', 'SGT', 'SF', 'SV', 'SIG', 'MT', 'PAC'
                ]

                variant_fields = [
                    variant.chromosome,
                    variant.position,
                    '.',
                    variant.ref_allele,
                    alt_alleles,
                    '.',
                    filter_fields,
                    info_fields,
                    ':'.join(format_fields),
                ]

                for sample in samples:

                    ploidy = variant.sample_ploidy[sample]

                    sample_fields = []

                    sample_fields.append(format_vcf_field(get_allele_values(
                        alleles, variant.strand_allele_counts[sample])))
                    sample_fields.append(format_vcf_field(
                        variant.allele_coverage_flags[sample]))
                    sample_fields.append(format_vcf_field(
                        variant.depth_flags[sample] or \
                        variant.filtered_sites_flags[sample]))
                    sample_fields.append(format_vcf_field(
                        variant.sample_ploidy[sample]))
                    sample_fields.append(get_vcf_genotype(
                        variant.sample_genotypes[sample], alleles, ploidy))
                    sample_fields.append(get_vcf_genotype(
                        variant.sample_subclonals[sample]['genotype'], alleles, ploidy))
                    sample_fields.append(format_vcf_field(
                        variant.sample_subclonals[sample]['frequency']))
                    sample_fields.append(format_vcf_field(
                        variant.sample_subclonal_valid_flag[sample]))

                    if sample != self.control_sample:
                        sample_fields.append(format_vcf_field(
                            variant.sample_significance_flags[sample]))
                        if variant.report_mutations[sample]:
                            sample_fields.append(get_vcf_mutations(
                                variant.sample_called_mutations[sample]))
                            sample_fields.append(format_vcf_field(
                                variant.sample_called_loh[sample]))
                        else:
                            sample_fields.extend(['.'] * 2)

                    variant_fields.append(':'.join(sample_fields))

                OUTPUT.write('\t'.join([format_value(x) for x in variant_fields]) + '\n')
