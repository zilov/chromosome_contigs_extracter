# Chromosome Contigs Extractor

Chromosome Contig Extractor is a Python script that extracts contigs aligned to chromosomes of a reference genome from an assembly file. It creates separate FASTA files for each chromosome, containing the aligned contigs. The script also identifies and reports chimeric contigs, i.e., contigs that align to multiple chromosomes.

## Features

- Extract aligned contigs from assembly file and save them in separate FASTA files for each chromosome
- Identify chimeric contigs and create a report with their headers, lengths, and the list of aligned chromosomes
- Customizable alignment length threshold
- Configurable output directory and file naming
- Logging for warnings and potential issues during the process

## Dependencies

- Python 3.6 or higher
- Pysam
- Click

## Installation

Clone the repository and install dependecies:

`bash
git clone https://github.com/yourusername/contig-extractor.git
pip3 install click && pip3 install pysam
`

## Usage

`bash
chromosome_contig_extractor.py assembly.fasta reference.fasta alignment.bam --alignment-threshold 100000 --output-prefix human --output-dir separate_chromosome_contigs --include-chimeric
`

This command will process the assembly and alignment files, creating separate FASTA files for each chromosome with the aligned contigs. Chimeric contigs will be included in the extracted contigs files and reported in a separate TSV file. The output files will be stored in the `separate_chromosome_contigs` directory with a prefix `human`.


## Arguments

- `<assembly_file>`: Path to the assembly file in FASTA format
- `<reference_file>`: Path to the reference genome file in FASTA format
- `<alignment_file>`: Path to the alignment file in SAM/BAM format
- `<threshold>`: (Optional) Minimum alignment length to consider (default: 50000)
- `<prefix>`: (Optional) Prefix for output file names (default: '')
- `<dir>`: (Optional) Directory to store output files (default: 'output')
- `--include-chimeric`: (Optional) Include chimeric contigs in extracted contigs files

## Outputs

- `<(prefix)_(chr_name)_aligned_contigs.fasta>` - list of fasta files for each reference genome chromosome
- `<(prefix)_possible_chimeric_contigs.tsv>` - tsv file which collects chimeric contigs headers, their lengths and list of chromosomes aligned

