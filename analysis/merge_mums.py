import os,sys
import numpy as np
import argparse

def parse_arguments():    
    parser = argparse.ArgumentParser(description="Plots a synteny plot of MUMs from mumemto")
    parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto', required=True)
    parser.add_argument('--gap','-g', dest='gap', help='maximum gap between MUMs (default: 5)', default=5, type=int)
    parser.add_argument('--num-mismatch','-x', dest='mismatch', help='maximum number of mismatch columns between MUMs', default=0, type=int)    
    args = parser.parse_args()
    
    if args.mismatch > args.gap:
        args.mismatch = 0
    assert args.gap > 0, "Gap must be greater than 0"
    return args

def merge_mums(args, mums):
    potential_glue = [mums[0]]
    blocks = []
    for idx in range(1, len(mums)):
        if mums[idx][2] != mums[idx - 1][2]:
            blocks.append(potential_glue)
            potential_glue = [mums[idx]]
            continue
        dist = (mums[idx][1] - mums[idx - 1][1] - mums[idx - 1][0])
        if dist[0] < args.gap and np.all(dist == dist[0]):
            potential_glue.append(mums[idx])
        else:
            blocks.append(potential_glue)
            potential_glue = [mums[idx]]
    return blocks

def write(args, merged):
    with open(os.path.splitext(args.mumfile)[0] + '_merged' + os.path.splitext(args.mumfile)[1], 'w') as f:
        for block in merged:
            blocklen = (block[-1][1][0] + block[-1][0]) - block[0][1][0]
            line = str(blocklen)
            starts = []
            for idx, strand in enumerate(block[0][2]):
                if strand == '+':
                    starts.append(block[0][1][idx])
                else:
                    starts.append(block[-1][1][idx])
            line += '\t' + ','.join([str(v) for v in starts]) + '\t' + ','.join(block[0][2])
            f.write(line + '\n')
    
def parse_mums(mumfile):
    mums = [l.split() for l in open(mumfile, 'r').read().splitlines()]
    mums = [(int(l[0]), np.array([int(v) for v in l[1].split(',')]), tuple(l[2].split(','))) for l in mums]
    mums = sorted(mums, key=lambda x: x[1][0])
    return mums

def main(args):
    mums = parse_mums(args.mumfile)
    merged = merge_mums(args, mums)
    avg_blocklen = np.mean([(block[-1][1][0] + block[-1][0]) - block[0][1][0] for block in merged])
    avg_mumlen = np.mean([m[0] for m in mums])
    nontrivial_blocks = sum([len(b) > 1 for b in merged])
    print('Average MUM length: %d' % avg_mumlen)
    print('Average merged MUM length: %d' % avg_blocklen)
    print('Average block len: %.3f' % np.mean([len(b) for b in merged]))
    print('Merged %d MUMs (%.2f%%)'%(nontrivial_blocks, nontrivial_blocks / len(mums) * 100))
    write(args, merged)
    
if __name__ == "__main__":
    args = parse_arguments()
    main(args)