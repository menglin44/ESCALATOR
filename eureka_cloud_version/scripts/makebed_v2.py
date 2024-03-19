import sys

infile=sys.argv[-2] # snpID, chr, pos, risk a, ref a, weight
outfile = open(sys.argv[-1], 'w')


for i, line in enumerate(open(infile)):
    line = line.strip().split()
    if not line[1].startswith('chr'):
        chr = 'chr' + line[1]
    else:
        chr = line[1]
    pos = line[2]
    if pos=='NA':
        continue
    pos1 = str(int(line[2])-1)
    pos2 = line[2]
    allele = line[3]
    ref = line[4]
    effect = line[5]
    snp = line[0]
    outfile.write('%s\n' % '\t'.join([chr, pos1, pos2, snp, allele, ref, effect]))

outfile.close()
