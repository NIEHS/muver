# -*- coding: utf-8 -*-

"""Console script for muver."""

import click

from allelic_fraction import get_allelic_fractions
from bias_distribution import calculate_bias_distribution_bam
from call_variants import call_variants as _call_variants
from depth_correction import write_corrected_bedgraph
from depth_distribution import (calculate_depth_distribution_bedgraph,
                                filter_regions_by_depth_bedgraph)
from pipeline import run_pipeline as _run_pipeline
from repeat_indels import fit_repeat_indel_rates
from repeats import create_repeat_file as _create_repeat_file
from utils import read_repeats
from wrappers.samtools import get_mpileup_depths


@click.group()
def main(args=None):
    pass

@main.command()
@click.option('--processes', '-p', default=1, type=int)
@click.option('--excluded_regions', default=None, type=str)
@click.argument('reference_assembly')
@click.argument('fastq_list')
@click.argument('control_sample')
@click.argument('experiment_directory')
def run_pipeline(reference_assembly, fastq_list, control_sample,
                 experiment_directory, processes, excluded_regions):
    _run_pipeline(
        reference_assembly,
        fastq_list,
        control_sample,
        experiment_directory,
        p=processes,
        excluded_regions=excluded_regions,
    )

@main.command()
@click.option('--chrom_sizes', default=None)
@click.option('--excluded_regions', default=None)
@click.argument('reference_assembly')
@click.argument('control_sample')
@click.argument('sample_list')
@click.argument('input_vcf')
@click.argument('output_header')
def call_variants(reference_assembly, control_sample, sample_list, input_vcf,
                  output_header, chrom_sizes, excluded_regions):
    _call_variants(
        reference_assembly,
        control_sample,
        sample_list,
        input_vcf,
        output_header,
        chrom_sizes=chrom_sizes,
        excluded_regions=excluded_regions,
    )

@main.command()
@click.argument('bam_file')
@click.argument('reference_assembly')
@click.argument('output_file')
def plot_allelic_fraction(bam_file, reference_assembly, output_file):
    get_allelic_fractions(
        bam_file,
        reference_assembly,
        output_file,
    )

@main.command()
@click.argument('fasta_file')
@click.argument('output_repeat_file')
def create_repeat_file(fasta_file, output_repeat_file):
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
@click.argument('chrom_sizes_file')
@click.argument('input_bedgraph')
@click.argument('output_bedgraph')
def depth_correction(y_int, scalar, mean_log, sd_log, slope, chrom_sizes_file,
                     input_bedgraph, output_bedgraph):
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
@click.argument('bam_file')
@click.argument('reference_assembly')
@click.argument('output_bedgraph')
def calculate_read_depths(bam_file, reference_assembly, output_bedgraph):
    get_mpileup_depths(
        bam_file,
        reference_assembly,
        output_bedgraph,
    )

@main.command()
@click.argument('bam_file')
@click.argument('reference_assembly')
@click.argument('output_bias_distribution')
def calculate_bias_distribution(bam_file, reference_assembly,
                                output_bias_distribution):
    calculate_bias_distribution_bam(
        bam_file,
        reference_assembly,
        output_bias_distribution,
    )

@main.command()
@click.option('--output_filtered_regions', type=str)
@click.argument('bedgraph_file')
@click.argument('chrom_sizes_file')
@click.argument('output_depth_distribution')
def calculate_depth_distribution(bedgraph_file, chrom_sizes_file,
                                 output_depth_distribution,
                                 output_filtered_regions):
    mu, sigma = calculate_depth_distribution_bedgraph(
        bedgraph_file,
        output_depth_distribution,
    )
    if output_filtered_regions:
        filter_regions_by_depth_bedgraph(
            bedgraph_file,
            chrom_sizes_file,
            mu,
            sigma,
            output_filtered_regions,
        )

@main.command()
@click.option('--output_plot_header', type=str)
@click.argument('bam_file')
@click.argument('repeats_file')
@click.argument('output_fits_file')
def calculate_repeat_indel_rates(bam_file, repeats_file, output_fits_file,
                                 output_plot_header):
    repeats = read_repeats(repeats_file)
    fit_repeat_indel_rates(
        repeats,
        bam_file,
        output_fits_file,
        output_plot_header,
    )

if __name__ == "__main__":
    main()
