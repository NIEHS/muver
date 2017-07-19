# -*- coding: utf-8 -*-

"""Console script for muver."""

import click

from allelic_fraction import get_allelic_fractions
from bias_distribution import calculate_bias_distribution_bam
from call_mutations import call_mutations as _call_mutations
from depth_correction import write_corrected_bedgraph
from depth_distribution import (calculate_depth_distribution_bedgraph,
                                filter_regions_by_depth_bedgraph)
from pipeline import run_pipeline as _run_pipeline
from reference import read_chrom_sizes_from_file
from repeat_indels import fit_repeat_indel_rates as _fit_repeat_indel_rates
from repeats import create_repeat_file as _create_repeat_file
from utils import read_repeats
from wrappers.samtools import get_mpileup_depths


@click.group()
def main(args=None):
    pass

@main.command()
@click.option('--processes', '-p', default=1, type=int,
              help='Number of processes to use.')
@click.option('--excluded_regions', default=None, type=click.Path(exists=True),
              help='Regions to exclude in mutation calling (BED format).')
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('fastq_list', type=click.Path(exists=True))
@click.argument('control_sample_name', type=str)
@click.argument('experiment_directory', type=str)
def run_pipeline(reference_assembly, fastq_list, control_sample_name,
                 experiment_directory, processes, excluded_regions):
    '''
    Run MuVer pipeline, starting with FASTQ files.

    FASTQ files and other parameters for each sample are specified in the
    FASTQ_LIST, a TXT file in tab-delimited format. For file specification,
    see tool documentation.

    Accepted fields for the FASTQ_LIST include:

    \b
    "Sample Name"
    "Mate 1 FASTQ"
    "Mate 2 FASTQ"
    "Ploidy"
    "CNV BedGraph"
    '''
    _run_pipeline(
        reference_assembly,
        fastq_list,
        control_sample_name,
        experiment_directory,
        p=processes,
        excluded_regions=excluded_regions,
    )

@main.command()
@click.option('--chrom_sizes', default=None)
@click.option('--excluded_regions', default=None, type=click.Path(exists=True),
              help='Regions to exclude in mutation calling (BED format).')
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('control_sample_name', type=str)
@click.argument('sample_list', type=click.Path(exists=True))
@click.argument('input_vcf', type=click.Path(exists=True))
@click.argument('output_header', type=str)
def call_mutations(reference_assembly, control_sample_name, sample_list,
                  input_vcf, output_header, chrom_sizes, excluded_regions):
    '''
    Call mutations from a HaplotypeCaller VCF file.

    Parameters for each sample are specified in the SAMPLE_LIST, a TXT file in
    tab-delimited format. For file specification, see tool documentation.

    Accepted fields for the SAMPLE_LIST include:

    \b
    "Sample Name"
    "Ploidy"
    "CNV BedGraph"
    "Strand Bias STD"
    "Filtered Sites"
    "Repeat Indel Fits"
    '''
    _call_mutations(
        reference_assembly,
        control_sample_name,
        sample_list,
        input_vcf,
        output_header,
        chrom_sizes=chrom_sizes,
        excluded_regions=excluded_regions,
    )

@main.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('output_file', type=str)
def plot_allelic_fraction(bam_file, reference_assembly, output_file):
    '''
    Plot allelic fractions over a BAM file.

    The OUTPUT_FILE contains a histogram of allelic fractions (tab-delimited
    TXT format).
    '''
    get_allelic_fractions(
        bam_file,
        reference_assembly,
        output_file,
    )

@main.command()
@click.argument('fasta_file', type=click.Path(exists=True))
@click.argument('output_repeat_file', type=str)
def create_repeat_file(fasta_file, output_repeat_file):
    '''
    Create repeat file for the FASTA_FILE sequence.
    '''
    _create_repeat_file(
        fasta_file,
        output_repeat_file,
    )

@main.command()
@click.argument('y_int', type=float)
@click.argument('scalar', type=float)
@click.argument('mean_log', type=float)
@click.argument('sd_log', type=float)
@click.argument('slope', type=float)
@click.argument('chrom_sizes_file', type=click.Path(exists=True))
@click.argument('input_bedgraph', type=click.Path(exists=True))
@click.argument('output_bedgraph', type=str)
def correct_depths(y_int, scalar, mean_log, sd_log, slope, chrom_sizes_file,
                     input_bedgraph, output_bedgraph):
    '''
    Correct values in a depth bedGraph file.

    Correction is performed using the sum of a log-normal cumulative
    distribution function and linear function.
    '''
    write_corrected_bedgraph(
        input_bedgraph,
        chrom_sizes_file,
        output_bedgraph,
        y_int,
        scalar,
        mean_log,
        sd_log,
        slope,
    )

@main.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('output_bedgraph', type=str)
def calculate_read_depths(bam_file, reference_assembly, output_bedgraph):
    '''
    Give depths in a BAM_FILE in bedGraph format.

    Depth values are calculated using the mpileup function of samtools.
    '''
    get_mpileup_depths(
        bam_file,
        reference_assembly,
        output_bedgraph,
    )

@main.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('output_bias_distribution', type=str)
def calculate_bias_distribution(bam_file, reference_assembly,
                                output_bias_distribution):
    '''
    Calculate distribution of strand biases.
    '''
    calculate_bias_distribution_bam(
        bam_file,
        reference_assembly,
        output_bias_distribution,
    )

@main.command()
@click.option('--output_filtered_regions', type=str, default=None,
              help='If OUTPUT_FILTERED_REGIONS is specified, positions to be '
              'filtered based on abnormal depth will be written to this file.')
@click.argument('bedgraph_file', type=click.Path(exists=True))
@click.argument('chrom_sizes_file', type=click.Path(exists=True))
@click.argument('output_depth_distribution', type=str)
def calculate_depth_distribution(bedgraph_file, chrom_sizes_file,
                                 output_depth_distribution,
                                 output_filtered_regions):
    '''
    Calculate distribution of depths in a bedGraph file.
    '''
    mu, sigma = calculate_depth_distribution_bedgraph(
        bedgraph_file,
        output_depth_distribution,
    )
    chrom_sizes = read_chrom_sizes_from_file(chrom_sizes_file)
    if output_filtered_regions:
        filter_regions_by_depth_bedgraph(
            bedgraph_file,
            chrom_sizes,
            mu,
            sigma,
            output_filtered_regions,
        )

@main.command()
@click.option('--output_plot_header', type=str, default=None,
              help='If OUTPUT_PLOT_HEADER is specified, PNGs are written '
              'that show the fit relative to the observed repeat indel rates.')
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('repeats_file', type=click.Path(exists=True))
@click.argument('output_fits_file', type=str)
def fit_repeat_indel_rates(bam_file, repeats_file, output_fits_file,
                                 output_plot_header):
    '''
    Fit a logistic function to log-transformed repeat indel rates.

    The OUTPUT_FITS_FILE gives the parameters for the derived fits.
    '''
    repeats = read_repeats(repeats_file)
    _fit_repeat_indel_rates(
        repeats,
        bam_file,
        output_fits_file,
        output_plot_header,
    )

if __name__ == "__main__":
    main()
