from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
import os
import argparse

def parse_arguments():    
    parser = argparse.ArgumentParser(description="Plots a synteny plot of MUMs from mumemto")
    # parser.add_argument('--filelist', '-f', dest='filelist', help='path to filelist from mumemto', required=True)
    # parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto', required=True)
    # parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist', required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input-prefix', '-i', dest='prefix', help='prefix for filelist, mums, and lengths files')
    group.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist')
    
    parser.add_argument('--filelist', '-f', dest='filelist', help='if the filelist is provided, then FASTA filenames are used as labels')
    parser.add_argument('--len-filter','-L', dest='lenfilter', help='only plot MUMs longer than threshold', default=0, type=int)
    parser.add_argument('--subsample','-s', dest='subsample', help='subsample every Nth mum', default=1, type=int)
    parser.add_argument('--center','-c', dest='center', action='store_true', help='center plot', default=False)
    parser.add_argument('--inversion-color','-ic', dest='inv_color', help='color for inversions', default='green')
    parser.add_argument('--mum-color','-mc', dest='mum_color', help='opacity of mums [0-1]', default=0.9, type=float)
    parser.add_argument('--fout','-o', dest='filename', help='plot fname (default: input_prefix)')
    parser.add_argument('--dims', dest='size', help='fig dimensions (inches) (default: 20,10)', default=(20,10), type=float, nargs=2)
    parser.add_argument('--dpi','-d', dest='dpi', help='dpi', default=500, type=int)
    
    args = parser.parse_args()
    if args.mumfile:
        args.prefix = os.path.splitext(args.mumfile)[0]
        args.lens = args.prefix + '.lengths'
    else:
        args.mumfile = args.prefix + '.mums'
        args.lens = args.prefix + '.lengths'
        
    if not args.filename:
        args.filename = args.prefix
    return args

def get_polygons(args, mums, genome_lengths, lenfilter=0, offset=0, color=(0.8, 0.8, 0.8), inv_color='green'):
    def points_to_poly(points):
        starts, ends = tuple(zip(*points))
        points = starts + ends[::-1]
        return points
    max_length = max(genome_lengths)
    centering = [0] * len(genome_lengths)
    if args.center:
        centering = [((max_length - g) / 2) + offset for g in genome_lengths]
    polygons = []
    colors = []    
    for (l, starts, strands) in mums:
        if l < lenfilter:
            continue
        inverted = strands[0] == '-'
        points = []
        for idx, (x, strand) in enumerate(zip(starts, strands)):
            points.append(((centering[idx] + x, idx), (centering[idx] + x + l, idx)) if strand == '+' else ((centering[idx] + genome_lengths[idx] - x - l, idx), (centering[idx] + genome_lengths[idx] - x, idx)))
            if not inverted and strand == '-':
                inverted = True
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
            elif inverted and strand == '+':
                inverted = False
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
        if len(points) > 2:
            polygons.append(points_to_poly(points))
            colors.append(color)
    return polygons, colors

def draw_synteny(genome_lengths, mums, lenfilter=0, dpi=500, size=None, genomes=None, filename=None):
    fig, ax = plt.subplots()
    max_length = max(genome_lengths)
    # centering = [0] * len(genome_lengths)
    # if args.center:
    #     centering = [(max_length - g) / 2 for g in genome_lengths]
    #     # ax.plot([centering[idx] + 0, centering[idx] + g], [idx, idx], color='gray', alpha=0.2, linewidth=0.5)
    # polygons = []
    # for (l, starts, strands) in mums:
    #     if l < lenfilter:
    #         continue
    #     points = [((centering[idx] + x, idx), (centering[idx] + x + l, idx)) if strand == '+' else ((centering[idx] + genome_lengths[idx] - x - l, idx), (centering[idx] + genome_lengths[idx] - x, idx)) for idx, (x, strand) in enumerate(zip(starts, strands))]
    #     starts, ends = tuple(zip(*points))
    #     points = starts + ends[::-1]
    #     polygons.append(points)
    polygons, colors = get_polygons(args, mums, genome_lengths, lenfilter=lenfilter, color=tuple([args.mum_color] * 3), inv_color=args.inv_color)
    ax.add_collection(PolyCollection(polygons, linewidths=.1, alpha=0.05, edgecolors=colors, facecolors=colors))
    ax.yaxis.set_ticks(range(len(genome_lengths)))
    ax.tick_params(axis='y', which='both',length=0)
    if genomes:
        ax.set_yticklabels(genomes)
    else:
        ax.yaxis.set_ticklabels([])
    ax.set_xlabel('bp')
    ax.set_ylabel('sequences')
    ax.set_ylim(0, len(genome_lengths))
    ax.set_xlim(0, max_length)
    fig.set_tight_layout(True)
    ax.axis('off')
    if size:
        fig.set_size_inches(*size)
    if filename:
        fig.savefig(os.path.join(os.path.dirname(args.mumfile), filename + ('' if filename.endswith('.png') else '.png')), dpi=dpi)
    return ax

def main(args):
    seq_lengths = [int(l.split()[1]) for l in open(args.lens, 'r').read().splitlines()]
    if args.filelist:
        genome_names = [os.path.splitext(os.path.basename(l.split()[0]))[0] for l in open(args.filelist, 'r').read().splitlines()]
    else:
        genome_names = None
    mums = parse_mums(args)
    draw_synteny(seq_lengths, mums, lenfilter=args.lenfilter, genomes=genome_names, filename=args.filename, dpi=args.dpi, size=args.size)

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