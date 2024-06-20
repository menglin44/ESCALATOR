import sys
import gzip

filename = sys.argv[-1]

if filename.endswith('.gz'):
    infile = gzip.open(filename)
else:
    infile = file(filename)


#trait='unknown'
pgs='unknown'
id_header='PGS ID'
for line in infile:
    if not line.startswith('#'):
        break
    #if 'Reported Trait' in line:
    #    trait = line.strip().split('=')[1].strip().replace(' ', '_')
    if 'format_version=2.0' in line:
        #print 'Header indicates the input file format as v2.0 from PGS catalog.'
        id_header = 'pgs_id'
    if id_header in line:
        pgs = line.strip().split('=')[1].strip().replace(' ', '_')
        
    
#assert 'unknown' not in [trait, pgs], "Didn't find [# Reported Trait] or [# PGS ID] in the header. Check your input file."
assert pgs!="unknown", "Didn't detect [# PGS ID] or [#pgs_id] in the header of the input file. Consider manually inputting a trait name."


print '%s' % (pgs)
