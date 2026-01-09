# -*- coding: utf-8 -*-
"""
this is the code used exclusivly for the logo plots in the publication
repeated 9925 to ensure the data sets they r drawn from are exactly the same in case
reviewers want all the supplementary data there is available.


"""

from Bio import AlignIO
import pandas as pd
import seqlogo
import logomaker as lm

#this doesnt need aligning because they are the exact same regex pull XX then CXC then XX
XXCXCXX_alignment = AlignIO.read(open(r"C:\PATH\XXCXCXX_9925.fasta"),'fasta')

#making and alignmen df for the XXCXCXX logo plot

def alnSiteCompositionDF(aln, characters="ACDEFGHIKLMNPQRSTVWY"):
    alnRows = aln.get_alignment_length()
    compDict = {char:[0]*alnRows for char in characters}
    for record in aln:
        header = record.id
        seq = record.seq
        for aaPos in range(len(seq)):
          aa = seq[aaPos]
          if aa in characters:
            compDict[aa][aaPos] += 1    
    return pd.DataFrame.from_dict(compDict)

# checking it did something after running the func

CXC_alignment_df = alnSiteCompositionDF(XXCXCXX_alignment)
CXC_alignment_df

#now I need to re index and normalize to 1

CXC_alignmentFREAK_df = CXC_alignment_df.div(CXC_alignment_df.sum(axis=1), axis=0)
CXC_alignmentFREAK_df
CXC_alignmentFREAK_df = CXC_alignmentFREAK_df.reset_index(drop=True)

#making the plot its not great for big stuff  tbh so you have to slice it:
    
lm.Logo(CXC_alignmentFREAK_df)

'''
thats it for the XXCXCXX next is the 4k alignment which i will have to drop indels for so the indel criteria for dropping is going to be 50% seqs have that position
ie 2000

restart the kernal between the above and below. thats how I ran it to ensure there was no carry over
'''
from Bio import AlignIO
import pandas as pd
import seqlogo
import logomaker as lm


alignment = AlignIO.read(open(r"C:\PATH\MauE_4k_unaligned_1cxc_3marks_9925_aligned.fasta"),'fasta')

# check the length of the alignment before removing indels
print("Alignment length %i" % alignment.get_alignment_length())
for record in alignment:
    print(record.seq + " " + record.id)
    
#making and alignmen df for the logo plot

def alnSiteCompositionDF(aln, characters="ACDEFGHIKLMNPQRSTVWY"):
    alnRows = aln.get_alignment_length()
    compDict = {char:[0]*alnRows for char in characters}
    for record in aln:
        header = record.id
        seq = record.seq
        for aaPos in range(len(seq)):
          aa = seq[aaPos]
          if aa in characters:
            compDict[aa][aaPos] += 1    
    return pd.DataFrame.from_dict(compDict)

# checking it did something after running the func

alignment_df = alnSiteCompositionDF(alignment)
alignment_df

columns = list(alignment_df.keys())
positions_dropped=[]
for ind in alignment_df.index:
    cut_off = 2000 #cut off im using for random 4k mauE
    sum_entries = 0
    for i in columns:
        sum_entries += alignment_df.at[ind,i]
    print(sum_entries)
    if sum_entries < cut_off:
        positions_dropped.append(ind)
        alignment_df = alignment_df.drop(ind)
        print('dropped '+ str(ind))

#now I need to re index and normalize to 1

alignmentFREAK_df = alignment_df.div(alignment_df.sum(axis=1), axis=0)
alignmentFREAK_df
alignmentFREAK_df = alignmentFREAK_df.reset_index(drop=True)
print(alignmentFREAK_df.shape)

start = 0
end = 26
for i in range(0,8):
    alignment_slice = alignmentFREAK_df.iloc[[i for i in range(start,end)]]
    lm.Logo(alignment_slice,figsize=(20,5))
    start += 26
    end += 26