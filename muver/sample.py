import csv
import glob
import os
from tempfile import NamedTemporaryFile

from repeat_indels import read_fits


class Sample(object):
    '''
    Object to hold information for a given sample in MuVer, including
    intermediate files.
    '''
    def __init__(self, sample_name, fastqs=[], kwargs=dict(),
                 exp_dir=None):

        self.sample_name = sample_name
        self.exp_dir = exp_dir
        self.fastqs = fastqs

        self.ploidy = 2

        for attr in (
            'cnv_bedgraph',
            'cnv_regions',
            'strand_bias_std',
            'filtered_sites',
            'strand_bias_distribution',
            'depth_distribution',
            'repeat_indel_header',
            'repeat_indel_fits',
            'repeat_indel_fits_dict',
            'merged_bam',
            'depth_bedgraph'
            '_mpileup_out',
        ):
            setattr(self, attr, None)

        for attr in (
            'realignment_logs',
            '_sams',
            '_same_chr_sams',
            '_mapq_filtered_sams',
            '_read_group_bams',
            '_deduplicated_bams',
            '_deduplication_metrics',
            '_interval_files'
            '_realigned_bams',
            '_fixed_mates_bams',
        ):
            setattr(self, attr, [])

        # SET SAMPLE ATTRIBUTES BASED ON KWARGS
        for field, value in kwargs.items():
            field = field.replace(' ', '_').lower().strip()
            if hasattr(self, field):

                if field == 'ploidy':
                    value = int(value.strip())
                if field == 'strand_bias_std':
                    value = float(value.strip())

                setattr(self, field, value)
            else:
                raise ValueError('{} is not an accepted field.'.format(field))

        if self.cnv_bedgraph:
            self.cnv_regions = self.read_cnv_bedgraph()

        if self.repeat_indel_fits:
            self.repeat_indel_fits_dict = read_fits(self.repeat_indel_fits)

    def read_cnv_bedgraph(self):
        '''
        Read in a bedGraph file to find ploidy at individual positions,
        used later to overwrite sample-wide ploidy.
        '''
        cnv_regions = dict()

        with open(self.cnv_bedgraph) as f:

            for line in f:

                chromosome, start, end, ploidy = line.strip().split()
                start = int(start) + 1
                end = int(end)
                ploidy = int(ploidy)

                for i in range(start, end + 1):
                    if ploidy != self.ploidy:
                        cnv_regions[(chromosome, i)] = ploidy

        return cnv_regions

    def generate_intermediate_files(self):
        '''
        Generate named temporary files, and generate names for output files
        based on the sample name. Both are then attributed to the object.
        '''
        # GENERATE INTERMEDIATE FILES FOR EACH FASTQ
        def named_temp(suffix=''):
            return NamedTemporaryFile(
                mode='w', dir=self.exp_dir, suffix=suffix)

        def named_file(file_dir, file_name):
            file_name = file_name.format(self.sample_name)
            return os.path.join(*(self.exp_dir, file_dir, file_name))

        if self.fastqs:
            n = range(len(self.fastqs))

            for attr, suffix in (
                ('_sams', '.sam'),
                ('_same_chr_sams', '.sam'),
                ('_mapq_filtered_sams', '.sam'),
                ('_read_group_bams', '.bam'),
                ('_deduplicated_bams', '.bam'),
                ('_deduplication_metrics', '.txt'),
                ('_realigned_bams', '.bam'),
                ('_fixed_mates_bams', '.bam'),
                ('_interval_files', '.intervals'),
            ):
                setattr(self, attr, [named_temp(suffix=suffix) for i in n])

            for attr, suffix in (
                ('_mpileup_out', '.txt'),
            ):
                setattr(self, attr, named_temp(suffix=suffix))

            realignment_logs = \
                ['{{}}.realignment.{}.log'.format(str(i)) for i in n]
            setattr(self, 'realignment_logs',
                    [named_file('logs', log) for log in realignment_logs])

            for attr, file_dir, file_name in (
                ('merged_bam', 'bams', '{}.bam'),
                ('depth_bedgraph', 'depth_distributions', '{}.bedGraph'),
                ('filtered_sites', 'filtered_sites',
                    '{}.filtered_sites.txt'),
                ('strand_bias_distribution', 'depth_distributions',
                    '{}.strand_bias_distribution.txt'),
                ('depth_distribution', 'depth_distributions',
                    '{}.depth_distribution.txt'),
                ('repeat_indel_fits', 'fits',
                    '{}.repeat_indel_fits.txt'),
                ('repeat_indel_header', 'fits',
                    '{}'),
            ):
                setattr(self, attr, named_file(file_dir, file_name))

    def get_intermediate_file_names(self):
        '''
        Iterate through object attributes and get temporary and output file
        names.  Return in a dict.
        '''
        file_names = dict()

        for attr, value in self.__dict__.items():

            if attr.startswith('_'):
                if isinstance(value, list):
                    file_names[attr] = [tmp.name for tmp in value]
                else:
                    file_names[attr] = value.name
            else:
                file_names[attr] = value

        return file_names

    def clear_temp_file_indices(self):
        '''
        Iterate through temporary files and delete any associated indices
        created by outside files.
        '''
        tempfile_names = []

        for attr, value in self.__dict__.items():
            if attr.startswith('_'):  # '_' identifies temporary files
                if isinstance(value, list):
                    tempfile_names.extend([tmp.name for tmp in value])
                else:
                    tempfile_names.append(value.name)

        for name in tempfile_names:
            for ext in ('.bai', '.idx'):
                to_remove = []

                to_remove.extend(
                    glob.glob('{}*{}'.format(name, ext)))
                to_remove.extend(
                    glob.glob('{}*{}'.format(os.path.splitext(name)[0], ext)))

                for fn in set(to_remove):
                    if fn != name:
                        os.remove(fn)


def read_samples_from_text(sample_info_fn, exp_dir=None):
    '''
    Iterate through an input file and use values to generate a list of sample
    objects.
    '''
    samples = []

    with open(sample_info_fn) as f:

        reader = csv.DictReader(f, delimiter='\t')
        fastqs = None

        for row in reader:

            if 'Mate 1 FASTQ' in row:
                fastqs = []
                if 'Mate 2 FASTQ' in row:
                    for f1, f2 in zip(row['Mate 1 FASTQ'].split(','),
                                      row['Mate 2 FASTQ'].split(',')):
                        fastqs.append((f1, f2))
                    row.pop('Mate 1 FASTQ')
                    row.pop('Mate 2 FASTQ')
                else:
                    for f1 in row['Mate 1 FASTQ'].split(','):
                        fastqs.append((f1))
                    row.pop('Mate 1 FASTQ')

            try:
                sample_name = row.pop('Sample Name').strip()
                samples.append(Sample(sample_name, fastqs=fastqs,
                               kwargs=row, exp_dir=exp_dir))
            except KeyError:
                raise ValueError('Sample Name required in sample file.')

    return samples


def generate_experiment_directory(exp_dir):
    '''
    For a given parent directory, create necessary subdirectories for storing
    output files if not found.
    '''
    for sub_dir in [
        'bams',
        'logs',
        'depth_distributions',
        'filtered_sites',
        'output',
        'gatk_output',
        'fits',
    ]:
        directory = os.path.join(exp_dir, sub_dir, '')

        if not os.path.isdir(directory):
            os.makedirs(directory)


def write_sample_info_file(samples, output_file):
    '''
    Write sample attributes to an output TXT file that may be read in later to
    generate a similar list of sample objects.
    '''
    sample_info = dict()

    fields = (
        'Sample Name',
        'Ploidy',
        'Filtered Sites',
        'Strand Bias Std',
        'Cnv Bedgraph',
        'Depth Distribution',
        'Strand Bias Distribution',
        'Merged Bam',
    )

    for sample in samples:
        sample_info[sample] = dict()
        for attr, value in sample.__dict__.items():
            field = attr.replace('_', ' ').title()
            if field in fields:
                if isinstance(value, list) or isinstance(value, tuple):
                    sample_info[sample][field] = ','.join(value)
                else:
                    sample_info[sample][field] = value

    with open(output_file, 'w') as OUT:
        writer = csv.DictWriter(OUT, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for sample in samples:
            writer.writerow(sample_info[sample])
