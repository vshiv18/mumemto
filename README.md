# **mumemto**: finding multi-MUMs and MEMs in pangenomes

<img src="img/polaroid_tattoo.png" alt="logo" width="292" align="left"/>

Mumemto identifies **maximal unique matches (multi-MUMs)** present across a collection of sequences. Multi-MUMs are defined as maximally matching substrings present in each sequence in a collection *exactly once*. Additionally, this tool can identify **multi-MEMs**, maximal exact matches present across sequences, without the uniqueness property. This method is uses the prefix-free parse (PFP) algorithm for suffix array construction on large, repetitive collections of text.

This tool uses PFP to efficiently identify multi-MUM/MEMs. Note that this applies only to highly repetitive texts (such as a collection of closely related genomes, likely intra-species such as a pangenome). We plan to support multi-MUM/MEM finding in more divergent sequences (inter-species, etc.) soon, however this would be less efficient without the PFP pre-processing step.

The base code from this repo was adapted from <a href="https://github.com/maxrossi91/pfp-thresholds">pfp-thresholds</a> repository written by <a href="https://github.com/maxrossi91">Massimiliano Rossi</a> and <a href="https://github.com/oma219/docprofiles">docprofiles</a> repository written by <a href="https://github.com/oma219">Omar Ahmed</a>. 

## Installation

### Docker/Singularity
Mumemto is available on `docker` and `singularity`:
```
### if using docker:
docker pull vshiv123/mumemto:latest
docker run vshiv123/mumemto:latest mumemto -h
### if using singularity:
singularity pull docker://vshiv123/mumemto:latest
./mumemto_latest.sif mumemto -h
```

### Compile from scratch
For starting out, use the commands below to download the repository and build the executable. After running the make command below,
the `mumemto` executable will be found in the `build/` folder. The following are dependencies: cmake, g++, gcc, libboost, zlib

```sh
git clone https://github.com/vshiv18/mumemto
cd mumemto

mkdir build 
cd build && cmake ..
make install
```

## Getting started

The basic workflow with `mumemto` is to compute the PFP over a collection of sequences, and identify multi-MUMs while computing the SA/LCP/BWT of the input collection. 

### Find multi-MUMs

```sh
mumemto mum -o <output_prefix> [input_fasta [...]]
```
Alternatively, you can find all multi-MEMs:
```sh
mumemto mem -o <output_prefix> [input_fasta [...]]
```

The command above takes in a list of fasta files as positional arguments and then generates output files using the output prefix. Alternatively, you can provide a file-list, which specifies a list of fastas and which document/class each file belongs in. Passing in fastas as positional arguments will auto-generate a filelist that defines the order of the sequences.

Use the `-h` flag to list the options for each mode: `mumemto mum -h`.
Mumemto mode options enable the computation of various different classes of exact matches:
<p align="center">
<img src="img/viz_def.png" alt="visual_guide" width="600" align="center"/>
</p>

`-k` allows for partial multi-MUM and MEMs (appearing in at least `N-k` sequences) and `--rare k` finds multi-MEMs that appear at most `k` times in each sequences (can be used with `-k` to find rare partial multi-MEMs).

**Format of the \*.mums file:**
```sh
[MUM length] [comma-delimited list of offsets within each sequence, in order of filelist] [comma-delimited strand indicators (one of +/-)]
```
The `*.mums` file contains each MUM as a separate line, where the first value is the match length, and the second is 
a comma-delimited list of positions where the match begins in each sequence. An empty entry indicates that the MUM was not found in that sequence (only applicable with *-k* flag). The MUMs are sorted in the output file
lexicographically based on the match sequence.

**Format of the \*.mems file:**
```sh
[MEM length] [comma-delimited list of offsets for each occurence] [comma-delimited list of sequence IDs, as defined in the filelist] [comma-delimited strand indicators (one of +/-)]
```
The `*.mems` file contains each MEM as a separate line with the following fields: (1) the match length, (2)
a comma-delimited list of offsets within a sequence, (3) the corresponding sequence ID for each offset given in (2). The MEMs are sorted in the output file
lexicographically based on the match sequence.

**Example of file-list file:**
```sh
/path/to/ecoli_1.fna 1
/path/to/salmonella_1.fna 2
/path/to/bacillus_1.fna 3
/path/to/staph_2.fna 4
```

## Visualization
<figure>
<img src="img/potato_syn_small.png" alt="potato_synteny"/>
<figcaption> <p align="center">Potato pangenome (assemblies from <a href='https://www.nature.com/articles/s41586-022-04822-x'>[Tang <i>et al.</i>, 2022]</a>)</p></figcaption>
</figure>
Mumemto can visualize multi-MUMs in a synteny-like format, highlighting conservation and genomic structural diversity within a collection of sequences.

After running `mumemto mum` on a collection of FASTAs, you can generate a visualization using:
```sh
/path/to/mumemto_repo/analysis/viz_mums.py (-i PREFIX | -m MUMFILE)
```
Use `viz_mums.py -h` to see options for customizability. As of now, only strict and partial multi-MUMs are supported (rare multi-MEM support coming soon).
