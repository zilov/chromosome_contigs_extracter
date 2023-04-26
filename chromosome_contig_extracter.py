#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#@created: 26.04.2023
#@author: Danil Zilov
#@contact: zilov.d@gmail.com

import os
import sys
import pysam
from collections import defaultdict
import click
import logging

@click.command()
@click.argument('assembly_file', type=click.Path(exists=True))
@click.argument('reference_file', type=click.Path(exists=True))
@click.argument('alignment_file', type=click.Path(exists=True))
@click.option('--alignment-threshold', default=50000, help='Minimum alignment length to consider a contig aligned to a chromosome.')
@click.option('--prefix', default='', help='A prefix to prepend to the output file names.')
@click.option('--outdir', default='./', type=click.Path(), help='The directory where output files will be saved.')
@click.option('--log-level', default='INFO', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']), help='The level of log messages to display.')
@click.option('--include-chimeric', is_flag=True, help='Include chimeric contigs in the extracted contigs files.')


def main(assembly_file, reference_file, alignment_file, alignment_threshold, prefix, outdir, log_level, include_chimeric):
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    extract_contigs(assembly_file, reference_file, alignment_file, alignment_threshold, prefix, outdir, include_chimeric)

def fasta_reader_yield(path_to_fasta_file):
    header = None
    with open(path_to_fasta_file) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header,"".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header,"".join(seq)

def write_fasta_file(output_file, contig_dict):
    with open(output_file, 'w') as fh:
        for header, sequence in contig_dict.items():
            fh.write(f'>{header}\n{sequence}\n')

def write_chimeric_contigs_file(output_file, chimeric_contigs, contig_sequences):
    with open(output_file, 'w') as f:
        f.write("contig_header\tcontig_length\tcsv_list_chromosome_names\n")
        for contig, chrom_set in chimeric_contigs.items():
            contig_length = len(contig_sequences[contig])
            f.write(f"{contig}\t{contig_length}\t{','.join(chrom_set)}\n")
            

def process_alignment_read(read, ref_chromosomes, contig_sequences, alignment_threshold, aligned_contigs, chimeric_contigs):
    if read.is_unmapped or read.reference_name not in ref_chromosomes or read.query_alignment_length < alignment_threshold:
        return
    contig_name = read.query_name
    if contig_name in contig_sequences:
        if contig_name in chimeric_contigs:
            chimeric_contigs[contig_name].add(read.reference_name)
        elif contig_name in aligned_contigs[read.reference_name]:
            return
        else:
            handle_new_alignment(contig_name, read.reference_name, aligned_contigs, chimeric_contigs)
    else:
        logging.warning(f"Contig {contig_name} not found in the assembly file. Skipping.")


def handle_new_alignment(contig_name, reference_name, aligned_contigs, chimeric_contigs):
    aligned_contigs[reference_name].add(contig_name)
    for chrom, contig_set in aligned_contigs.items():
        if chrom != reference_name and contig_name in contig_set:
            aligned_contigs[chrom].remove(contig_name)
            chimeric_contigs[contig_name].add(chrom)
            chimeric_contigs[contig_name].add(reference_name)
            logging.warning(f"Possible chimeric contig detected: {contig_name}, aligned to {','.join(chimeric_contigs[contig_name])}")
            break


def extract_contigs(assembly_file, reference_file, alignment_file, alignment_threshold, output_prefix, output_dir, include_chimeric):
    ref_chromosomes = dict(fasta_reader_yield(reference_file))
    contig_sequences = dict(fasta_reader_yield(assembly_file))

    aligned_contigs = defaultdict(set)
    chimeric_contigs = defaultdict(set)

    with pysam.AlignmentFile(alignment_file, 'rb') as af:
        for read in af.fetch(until_eof=True):
            process_alignment_read(read, ref_chromosomes, contig_sequences, alignment_threshold, aligned_contigs, chimeric_contigs)

    os.makedirs(output_dir, exist_ok=True)

    for chrom, contig_set in aligned_contigs.items():
        if include_chimeric:
            contig_set.update({contig for contig in chimeric_contigs if chrom in chimeric_contigs[contig]})
        output_file = os.path.join(output_dir, f"{output_prefix}{chrom}_aligned_contigs.fasta")
        contig_dict = {contig: contig_sequences[contig] for contig in contig_set}
        write_fasta_file(output_file, contig_dict)
        
    chimeric_contigs_file = os.path.join(output_dir, f"{output_prefix}_possible_chimeric_contigs.tsv")
    write_chimeric_contigs_file(chimeric_contigs_file, chimeric_contigs, contig_sequences)



if __name__ == '__main__':
    main()