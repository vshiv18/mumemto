import sys, os
import shutil
import subprocess
from tqdm.auto import tqdm
import argparse

def parse_arguments():    
    parser = argparse.ArgumentParser(description="Converts MUMs to parsnp format and runs the parsnp aligner")
    parser.add_argument('--filelist', '-f', dest='filelist', help='path to filelist from mumemto', required=True)
    parser.add_argument('--mums', '-m', dest='mums', help='path to *.mum file from mumemto', required=True)
    parser.add_argument('--lens','-l', dest='lens', help='lengths file, first column is seq length in order of filelist', required=True)
    parser.add_argument('--dir','-o', dest='wkdir', help='parsnp working directory', default='parsnp')
    parser.add_argument('--parsnp-path', dest='parsnp_path', help='parsnp exec path', default='~/software/parsnp_msa/parsnp')
    parser.add_argument('--threads', '-t', dest='threads', help='number of threads to run parsnp on', default=1)
    
    # parser.add_argument('--parsnp-path', dest='parsnp_path', help='parsnp exec path', default='~/software/parsnp_msa/parsnp')
    args = parser.parse_args()
    return args

def lcp_to_parsnp_coords(p, mumlen, strand, seqlen):
    return p if strand == '+' else seqlen - mumlen - p
def mum_lcp_to_parsnp(mum, lens):
    return (mum[0], tuple([lcp_to_parsnp_coords(p, mum[0], s, l) for (p, s, l) in zip(mum[1], mum[2], lens)]), mum[2])

inifile = """;Parsnp configuration File
;
[Reference]
file={0}
reverse=0
[Query]
{1}
[MUM]
anchors=20
anchorfile={3}/{2}
anchorsonly=0
calcmumi=0
mums=20
mumfile=
filter=1
factor=2.0
extendmums=0
[LCB]
recombfilter=0
cores=48
diagdiff=0.12
doalign=2
c=21
d=300
q=30
p=15000000
icr=0
unaligned=0
[Output]
outdir={3}
prefix=parsnp
showbps=1
"""

def main(args):
    os.makedirs(args.wkdir, exist_ok=True)
    
    files = [f.split()[0] for f in open(args.filelist, 'r').read().splitlines()]
    ref, files = files[0], files[1:]
    shutil.copyfile(ref, os.path.join(args.wkdir, os.path.basename(ref)+'.ref'))
    ref = os.path.join(args.wkdir, os.path.basename(ref)+'.ref')
    all_lens = [int(l.split()[1]) for l in open(args.lens, 'r').read().splitlines()]

    outfile = open(os.path.join(args.wkdir, 'parsnp_formatted.mums'), 'w')
    ### reformat mums
    for l in open(args.mums, 'r').readlines():
        l = l.strip().split()
        l = (int(l[0]), tuple([int(v) for v in l[1].split(',')]), tuple(l[2].split(',')))
        l = mum_lcp_to_parsnp(l, all_lens)
        outfile.write(f"{l[0]}\t{','.join(map(str, l[1]))}\t{','.join(l[2])}\n")
    outfile.close()
    # if not os.path.exists(args.mums.replace('.mums', '.fai')):
    #     os.system('samtools faidx ' + args.mums.replace('.mums', ''))

    inilines = "\n".join(['file{0}={1}\nreverse{0}=0'.format(x + 1, f) for x, f in enumerate(files)])

    with open(os.path.join(args.wkdir, 'parsnpAligner.ini'), 'w') as out:
        out.write(inifile.format(ref, inilines, 'parsnp_formatted.mums', args.wkdir))
    
    command = '/usr/bin/time -v {0} -r {1} -d {2} -c -o {3} -v -p {4} -i {3}/parsnpAligner.ini --skip-phylogeny'.format(args.parsnp_path, ref, ' '.join(files), args.wkdir, args.threads)
    print('Running:', command)
    subprocess.run(command, shell=True)
    
if __name__ == "__main__":
    args = parse_arguments()
    main(args)