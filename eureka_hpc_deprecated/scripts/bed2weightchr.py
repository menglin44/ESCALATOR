'''
This script match .bed format (weight file) against .pvar file, 
1) flip any loci if need to be (from bed). 
2) remove mismatched loci (from bed).
3) rename snpID in bed to match pvar

Usage:
python bed2weightchr.py [T/F] [BED file] [pvar file] [output file]
Outputs a weight file for PLINK: snpID, risk allele, weight

If long version chosen ("F"), then also outputs record files in semi-bed file format - chr, bp, snpID(as in intput bed), risk allele, weight
2) ***_flipped.list: bed variants that got flipped on the strand
3) ***_mismatch.list: bed variants whose allele codes dont match w/ pvar, and are removed
4) ***_missing.list: bed variants whose positions are not found in pvar, and are removed
5) output a more complete version of variants that ended up being used: chr, bp, OriginalSNPID, UpdatedSNPID, RiskA, RefA, Weight
'''


import sys
import time
import numpy as np

speed = sys.argv[-4]
bed = sys.argv[-3] # small file
pvar = sys.argv[-2] # big file
outfile = file(sys.argv[-1], 'w') #  weight file format for PLINK --score, as snpID, allele code, weight


starttime=time.time()

# helper dict used for flipping
flip = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

# helper function to flip strand for >=1 character strings
def flipstr(geno):
    flipped = []
    for code in geno:
        flipped.append(flip.get(code, 'n'))
    return ''.join(flipped)


# build a dictionary of positions:sorted_alleles from bed
# no worries on indel misnamed as snp vice versa - will be removed if allele codes not matched
# also build a dictionary of position_sortedalleles:BEDsnpID_riskAllele_refAllele_weight
beddict={}
posallele2weight={}
for i,line in enumerate(file(bed)):
    if i % 10000==0:
        print 'on bed line %s' % i
    line = line.strip().split() # chr, bp-1, bp, snpID, riskA, refA, weight
    pos = line[0].replace('chr', '') + ':' + line[2]
    alleles = '_'.join(sorted(line[4:6]))
    beddict[pos] = alleles
    snpID = line[3].replace('_', '$') #some line[3] already has _ in its name, creating a problem for downstream split and extract
    posallele2weight[pos+'_'+alleles] = snpID + '_' + line[4] + '_' + line[5] + '_' + line[6] 

bedpos = set(beddict.keys())
dicttime=time.time()
print 'Finished building the dictionary from bed file: %ss' % str(dicttime - starttime)

#preprare a list of position_sortedalleles in BED file that has matched positions with pvar (regardless of allele code)
usesnps=[]
# prepare a list of position_sortedalleles in BED files that has matched allele codes (flipped or not) with pvar
codesnps=[]

# read in pvar an compare
pastheader=False
if speed == 'T': # brief version, no output of variant list being filtered, flipped etc.
    print 'Brief version chosen. Wont output extra documents of variants flipped, unmatched, or missing.'
    for i,line in enumerate(file(pvar)):
        if i % 10000 == 0:
            print 'on pvar line %s' % i
        if pastheader:
            line = line.strip().split()
            pos = line[0] + ':' + line[1]
            alleles = '_'.join(sorted(line[3:5]))
            snp = line[2]
            if pos not in bedpos: #position not found in BED
                continue
            elif beddict.get(pos, 'n_n') == alleles: # positions overlap and allele code matches between bed and pvar
                tempinfo = posallele2weight.get(pos+'_'+alleles, 'NA_NA_NA_NA').split('_') # OriginalSnpID, riskA, refA, weight
                if tempinfo[1] == line[4]: #risk allele is also the alternative allele:
                    if tempinfo[3] =='NA':
                        tempinfo[3] = '0'
                    outfile.write('%s\n' % '\t'.join([snp, tempinfo[1], tempinfo[3]]))
                else: # risk allele is the reference allele in pvar
                    if tempinfo[3]=='NA':
                        weight='0'
                    elif tempinfo[3].startswith('-'): #negative value
                        weight = tempinfo[3].strip('-')
                    else: # positive value,
                        weight = '-' + tempinfo[3]
                    outfile.write('%s\n' % '\t'.join([snp, tempinfo[2], weight]))
            elif beddict.get(pos, 'n_n') == '_'.join(sorted([flipstr(line[3]), flipstr(line[4])])): # positions overlap but strand flipped
                bed_alleles = beddict.get(pos, 'n_n')
                tempinfo = posallele2weight.get(pos+'_'+ bed_alleles, 'NA_NA_NA_NA').split('_') # OriginalSnpID, riskA, refA, weight
                #flipped = flip.get(tempinfo[1], 'unknown_code')
                flipped = flipstr(tempinfo[1])
                if flipped == line[4]: #risk allele is also the alternative allele
                    if tempinfo[3] =='NA':
                        tempinfo[3] = '0'
                    outfile.write('%s\n' % '\t'.join([snp, flipped, tempinfo[3]]))
                else: # risk allele is the reference allele in pvar
                    if tempinfo[3]=='NA':
                        weight='0'
                    elif tempinfo[3].startswith('-'): #negative value
                        weight = tempinfo[3].strip('-')
                    else: # positive value
                        weight = '-' + tempinfo[3]
                    flipped = flipstr(tempinfo[2])
                    outfile.write('%s\n' % '\t'.join([snp, flipped, weight]))
            else: # position overlap, but allele code don't match regardless of flipping
                continue
        elif line.startswith('#CHROM'):
            pastheader=True    
else: # long version, output lists of variants being flipped on strand, or discarded due to either mismatching code or bp
    out_flip = file(bed.replace('.bed', '_flipped.list'), 'w') # flip and save in weight file
    out_mismatch = file(bed.replace('.bed', '_mismatch.list'), 'w') # need to remove from weight
    out_missing = file(bed.replace('.bed', '_missing_in_pvar.list'), 'w') # those in BED that are not found in pvar
    out_used = file(bed.replace('.bed', '_cleaned_forRecord.list'), 'w') # a record of variants being used
    for i,line in enumerate(file(pvar)):
        if i % 10000 ==0:
            print 'on pvar line %s' % i
        if pastheader:
            line = line.strip().split()
            pos = line[0] + ':' + line[1]
            alleles = '_'.join(sorted(line[3:5]))
            snp = line[2]
            if pos not in bedpos: #position not found in BED
                continue
            elif beddict.get(pos, 'NA') == alleles: # positions overlap and allele code matches between bed and pvar
                tempinfo = posallele2weight.get(pos+'_'+alleles, 'NA_NA_NA_NA').split('_') # OriginalSnpID, riskA, refA, weight
                if tempinfo[1] == line[4]: #risk allele is also the alternative allele:
                    if tempinfo[3] =='NA':
                        tempinfo[3] = '0'
                    outfile.write('%s\n' % '\t'.join([snp, tempinfo[1], tempinfo[3]]))
                    out_used.write('%s\n' % '\t'.join(pos.split(':') + [tempinfo[0].replace('$', '_'), snp] + tempinfo[1:]))
                else: # risk allele is the reference allele in pvar
                    if tempinfo[3] == 'NA':
                        weight = '0'
                    elif tempinfo[3].startswith('-'): #negative value
                        weight = tempinfo[3].strip('-')
                    else: # positive value
                        weight = '-' + tempinfo[3]
                    outfile.write('%s\n' % '\t'.join([snp, tempinfo[2], weight]))
                    out_used.write('%s\n' % '\t'.join(pos.split(':') + [tempinfo[0].replace('$', '_'), snp, tempinfo[2], tempinfo[1], weight]))
                usesnps.append(pos + '_' + alleles)
                codesnps.append(pos + '_' + alleles)  
            elif beddict.get(pos, 'n_n') == '_'.join(sorted([flipstr(line[3]), flipstr(line[4])])): # positions overlap but strand flipped
                bed_alleles = beddict.get(pos, 'n_n')
                tempinfo = posallele2weight.get(pos+'_'+bed_alleles, 'NA_NA_NA_NA').split('_') # OriginalSnpID, riskA, refA, weight
                flipped = flipstr(tempinfo[1])
                if flipped == line[4]: #risk allele is also the alternative allele
                    if tempinfo[3] =='NA':
                        tempinfo[3] = '0'
                    outfile.write('%s\n' % '\t'.join([snp, flipped, tempinfo[3]]))
                    out_used.write('%s\n' % '\t'.join(pos.split(':') + [tempinfo[0].replace('$', '_'), snp, flipped, flipstr(tempinfo[2]), tempinfo[3]])) # needed to flip both risk and ref for the correct record
                else: # risk allele is the reference allele in pvar
                    if tempinfo[3] == 'NA':
                        weight = '0'
                    elif tempinfo[3].startswith('-'): #negative value
                        weight = tempinfo[3].strip('-')
                    else: # positive value
                        weight = '-' + tempinfo[3]
                    flipped = flipstr(tempinfo[2])
                    outfile.write('%s\n' % '\t'.join([snp, flipped, weight]))
                    out_used.write('%s\n' % '\t'.join(pos.split(':') + [tempinfo[0].replace('$', '_'), snp, flipped, flipstr(tempinfo[1]), weight])) # needed to flip both risk and ref for the correct record
                out_flip.write('%s\n' % '\t'.join(pos.split(':') + [tempinfo[0].replace('$', '_'), snp] + tempinfo[1:])) # write chr, pos, OriginalsnpID, pvarsnpID, riskAllele(BED), refAllele(BED), weight for a record
                usesnps.append(pos+'_'+bed_alleles)
                codesnps.append(pos+'_'+bed_alleles)    
            else: # position overlap, but allele codes don't match regardless of flipping
                bed_alleles = beddict.get(pos, 'n_n')
                #tempinfo = posallele2weight.get(pos+'_'+bed_alleles, 'NA_NA_NA_NA').split('_') # OriginalSnpID, riskA, refA, weight
                #out_mismatch.write('%s\n' % '\t'.join(pos.split(':') + [tempinfo[0], snp] + tempinfo[1:])) # chr, bp, OriginalsnpID, pvarsnpID, riskAllele(BED), refAllele(BED), weight for a record
                usesnps.append(pos+'_'+bed_alleles)
        elif line.startswith('#CHROM'):
            pastheader=True    
    # output the list of BED positions unmatched in pvar, or mismatched allele codes in pvar
    bedvar = np.array(posallele2weight.keys())
    usedvar = np.array(usesnps)
    codevar = np.array(codesnps)
    missing_var = np.setdiff1d(bedvar, usedvar)
    mismatched_var = np.setdiff1d(usedvar, codevar)
    for var in missing_var:
        tempinfo = posallele2weight.get(var, 'NA_NA_NA_NA').split('_')
        out_missing.write('%s\n' % '\t'.join(var.split('_')[0].split(':') + [tempinfo[0].replace('$', '_')] + tempinfo[1:])) # chr, bp, OriginalsnpID, riskAllele(BED), refAllele(BED), weight for a record
    
    for var in mismatched_var:
        tempinfo = posallele2weight.get(var, 'NA_NA_NA_NA').split('_')
        out_mismatch.write('%s\n' % '\t'.join(var.split('_')[0].split(':') + [tempinfo[0].replace('$', '_')] + tempinfo[1:]))
        
    out_flip.close()
    out_mismatch.close()
    out_missing.close()
    out_used.close()

outfile.close()
pvartime=time.time()
print 'Finished reading and comparing pvar against dict: %ss' % str(pvartime - dicttime)




