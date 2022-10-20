#!/usr/bin/env python
"""

:Author: Ji Research Group/Stanford Genome Technology Center
:Contact: sgrimes@stanford.edu
:Creation date: 07/08/2021
:Description: 

This script determines genotype per read, at a supplied genome region.
Current usage is to assess CRISPR base editing
    
Revisions:
    mm/dd/yyyy: Description

"""
import sys, os, re, pysam, csv, numpy as np

script_name = os.path.basename(__file__)
print("Running ", script_name)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Files are assumed to be in current path (or fully qualified path is given)
if len(sys.argv) < 3:
  print("Usage: ", script_name, "<sam_file|bam_file> <region> <out_prefix>")
  sys.exit(1)

sam_file = sys.argv[1]
if os.path.isfile(sam_file) and os.access(sam_file, os.R_OK):
  sam_input = pysam.Samfile(sam_file,'rb') if sam_file[-3:] == 'bam' else pysam.Samfile(sam_file,'r')
  sam_fname = os.path.basename(sam_file)
else:
  print("Unable to open sam/bam file for input:", sam_file)
  sys.exit(1)
  
arg_coords = sys.argv[2].split(':')
target_chr = arg_coords[0]
target_coords = arg_coords[1].split('-')

if len(sys.argv) > 3:
  oprefix = sys.argv[3]
else:
  oprefix = sam_fname[0:-4]
  
try:
  out_file = oprefix + '.be_genotypes.txt'
  transcript_out = open(out_file, 'w')
except:
  print("Unable to transcript/genotype file for output: ", out_file)
  sys.exit(1)

CHKRD_NAMES = ['m54278_190912_162615/14549522/ccs', 'm54278_190912_162615/37094193/ccs']

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define internal modules                                                     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#   
def printDebug(i, qname, chkrdnames):
  return True if (i <1 or qname in chkrdnames) else False; 
  
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Read sam file and check if spans complete genomic coordinate range (from .bed)
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
i = 0; tot_rds = 0;
out_csv = csv.writer(transcript_out, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
#out_csv.writerow(['Qname', 'Genotype'])

target_start = int(target_coords[0])
target_end = int(target_coords[1])
ref_coords = [*range(target_start, target_end, 1)]

for bam_rd in sam_input.fetch(region=sys.argv[2]):
   i += 1
   rd_positions = bam_rd.get_aligned_pairs()
   target_positions = [pos_tuple for pos_tuple in rd_positions if pos_tuple[1] in ref_coords]
   target_bases = ''
   for coord in ref_coords:
     query_pos = [pos_tuple[0] for pos_tuple in target_positions if pos_tuple[1] == coord]
     if len(query_pos) == 0:
       query_base = 'N'
     elif query_pos[0] == None:
       query_base = '*'
     else:
       query_base = bam_rd.query_sequence[query_pos[0]]
     target_bases = target_bases + query_base
   out_csv.writerow([bam_rd.qname, target_bases])

   #if i > 250: break
   
sam_input.close()
transcript_out.close()
