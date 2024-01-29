# **Mumma Mia**: finding multi-MUMs using prefix-free parsing for big BWTs

This code is based on the [pfp-thresholds](https://github.com/maxrossi91/pfp-thresholds) repository written by [Massimiliano Rossi](https://github.com/maxrossi91) and [docprofiles] repository written by [Omar Ahmed](https://github.com/oma219). 

This repository identifies **maximal unique matches (multi-MUMs)** present across a collection of sequences. Multi-MUMs are defined as maximally matching substrings present in each sequence in a collection *exactly once*. This method is uses the prefix-free parse (PFP) algorithm for suffix array construction on large, repetitive collections of text.

This tool uses PFP to efficiently identify multi-MUMs. Note that this applies only to highly repetitive texts (such as a collection of closely related genomes, likely intra-species such as a pangenome). We plan to support multi-MUM finding in more divergent sequences (inter-species, etc.) soon, however this would be less efficient without the PFP pre-processing step.

## Installation

For starting out, use the commands below to download the repository and build the executable. After running the make command below,
the `pfp_mum` executable will be found in the `build/` folder.

```sh
git clone git@github.com:vshiv18/pfp-mum.git
cd pfp-mum

mkdir build 
cd build && cmake ..
make install

export PFPMUM_BUILD_DIR=$(pwd)
```

## Getting started

The basic workflow with `pfp_mum` is to compute the PFP over a collection of sequences, and identify multi-MUMs while computing the SA/LCP/BWT of the input collection. 

### Find multi-MUMs

```sh
pfp_mum -f <input_list> -o <output_prefix>
```

The command above takes in a file-list of multiple genomes and then generates output files using the output prefix. In the file-list, you can specify a list of 
genomes and then specify which document/class each genome belongs in.

**Example of file-list file:**
```sh
/path/to/ecoli_1.fna 1
/path/to/ecoli_2.fna 1
/path/to/salmonella_1.fna 2
/path/to/bacillus_1.fna 3
/path/to/staph_2.fna 4
```

**Format of the *.mums file:**
```sh
[MUM length] [comma-delimited list of offsets within each sequence, in order of filelist]
```
The `*.mums` file contains each MUM as a separate line, where the first value is the match length, and the second is 
a comma-delimited list of positions where the match begins in each sequence. The MUMs are sorted in the output file
lexicographically based on the match sequence.
