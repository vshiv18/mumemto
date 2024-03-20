from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
import os
import argparse

def parse_arguments():    
    parser = argparse.ArgumentParser(description="Plots a synteny plot of MUMs from mumemto")
    parser.add_argument('--filelist', '-f', dest='filelist', help='path to filelist from mumemto', required=True)
    parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto', required=True)
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist', required=True)
    parser.add_argument('--len-filter','-L', dest='lenfilter', help='only plot MUMs longer than threshold', default=0, type=int)
    parser.add_argument('--subsample','-s', dest='subsample', help='subsample every Nth mum', default=1, type=int)
    parser.add_argument('--center','-c', dest='center', action='store_true', help='center plot')
    parser.add_argument('--fout','-o', dest='filename', help='plot fname', default='mums')
    
    # parser.add_argument('--parsnp-path', dest='parsnp_path', help='parsnp exec path', default='~/software/parsnp_msa/parsnp')
    args = parser.parse_args()
    return args

def draw_synteny(genome_lengths, mums, lenfilter=0, dpi=500, size=None, genomes=None, filename=None):
    fig, ax = plt.subplots()
    max_length = max(genome_lengths)
    centering = [0] * len(genome_lengths)
    for idx, g in enumerate(genome_lengths):
        if args.center:
            centering[idx] = (max_length - g) / 2
        ax.plot([centering[idx] + 0, centering[idx] + g], [idx, idx], color='gray', alpha=0.2, linewidth=0.5)
    polygons = []
    for (l, starts, strands) in mums:
        if l < lenfilter:
            continue
        points = [((centering[idx] + x, idx), (centering[idx] + x + l, idx)) if strand == '+' else ((centering[idx] + x, idx), (centering[idx] + x - l, idx)) for idx, (x, strand) in enumerate(zip(starts, strands))]
        starts, ends = tuple(zip(*points))
        points = starts + ends[::-1]
        polygons.append(points)
    ax.add_collection(PolyCollection(polygons, linewidths=0))
    ax.yaxis.set_ticks(range(len(genome_lengths)))
    if genomes:
        ax.set_yticklabels(genomes)
    ax.set_xlabel('bp')
    fig.set_tight_layout(True)
    ax.set_ylabel('genomes')
    fig.set_dpi(dpi)
    ax.axis('off')
    if size:
        fig.set_size_inches(*size)
    else:
        fig.set_size_inches((6.4, len(genome_lengths) // 5))
    if filename:
        fig.savefig(os.path.join(os.path.dirname(args.mumfile), filename + ('' if filename.endswith('.png') else '.png')))
    return ax

def main(args):
    seq_lengths = [int(l.split()[1]) for l in open(args.lens, 'r').read().splitlines()]
    genome_names = [os.path.splitext(os.path.basename(l.split()[0]))[0] for l in open(args.filelist, 'r').read().splitlines()]
    mums = parse_mums(args)
    draw_synteny(seq_lengths, mums, lenfilter=args.lenfilter, genomes=genome_names, filename=args.filename)

def parse_mums(args):
    count = 0
    for l in open(args.mumfile, 'r').readlines():
        if count % args.subsample == 0:
            l = l.strip().split()
            yield int(l[0]), tuple([int(v) for v in l[1].split(',')]), tuple(l[2].split(','))
        count += 1
if __name__ == "__main__":
    args = parse_arguments()
    main(args)