import sys

bed = sys.argv[-2]
outfile = file(sys.argv[-1], 'w')

at_cg = file(bed.replace('.bed', '_at_cg.list'), 'w')


for line in file(bed):
    line = line.strip().split()
    if sorted(line[4:6]) in (['A','T'], ['C', 'G']):
        at_cg.write('%s\n' % '\t'.join(line))
    else:
        outfile.write('%s\n' % '\t'.join(line))

outfile.close()
at_cg.close()
