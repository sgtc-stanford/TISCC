#!/usr/bin/env python
import sys, os, re, pysam, csv, gzip, string, distance
import numpy as np, pandas as pd

DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA") 
A1_3PRIME = 'CTGGCTCCTCTGTATGTTGGAGAAT'
SI_5PRIME = 'AAGCAGTGGTATCAACGCAGAGTAC'
FWD_10X_ADAPTER = 'TTTCTTATATGGG'
REV_10X_ADAPTER = 'CCCATATAAGAAA' 
ILLUMINA_R1_FWD = 'CTACACGACGCTCTTCCGATCT'
ILLUMINA_R1_REV = 'AGATCGGAAGAGCGTCGTGTAG'

ADAPTER_3PRIME = [A1_3PRIME, A1_3PRIME[::-1], A1_3PRIME[::-1].translate(DNA_COMPLEMENT)] 
ADAPTER_5PRIME = [SI_5PRIME, SI_5PRIME[::-1], SI_5PRIME[::-1].translate(DNA_COMPLEMENT)] 

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define reusable modules to be imported/called from long read barcode-matching scripts     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+# 

def extract_softclips(samrd,offset=0):
  soft_clip_fwd = 'NNNNNNNNN'; soft_clip_rev = 'NNNNNNNNN';
  sam_cigar = samrd.cigartuples
  if sam_cigar and sam_cigar[0][0] == 4:
    soft_clip_len = sam_cigar[0][1]
    soft_clip_fwd = samrd.query_sequence[:soft_clip_len+offset] 
  if sam_cigar and sam_cigar[-1][0] == 4:
    soft_clip_len = sam_cigar[-1][1]
    soft_clip_rev = samrd.query_sequence[(0-soft_clip_len-offset):]
  return {'fwd': soft_clip_fwd, 'rev': soft_clip_rev}
  
def first_aligned_bases(samrd, n=10):
  qstart = samrd.query_alignment_start
  qend = samrd.query_alignment_end
  return {'fwd': samrd.query_sequence[qstart:qstart+n], 'rev': samrd.query_sequence[qend-n:qend]}

def reverse_complement(seq):
  return seq.translate(DNA_COMPLEMENT)[::-1]

def closest_pos_match(pos_list, pos):      
  pos_array = np.asarray(pos_list)
  idx = (np.abs(pos_array - pos)).argmin()
  return pos_array[idx] 
  
def calc_edit_distance(barcode, softclip_seq, search_start):
  dist = []; bc_len=len(barcode); sc_len=len(softclip_seq);
  
  #First look for exact match
  #To allow for multiple matches: 
  #  all_matches = [m.start(0) for m in re.finditer(barcode, search_seq)] 
  #  if no matches: all_matches list will be empty
  match_pos = softclip_seq.find(barcode)
  if match_pos > -1:
    return [0, match_pos]
  
  #If no exact match, search starting at expected position, radiating outwards
  # (as defined in search_pos_order)
  else:
    i_range = search_pos_order(sc_len, search_start, bc_len)
    for i in i_range:
      dist.append(distance.levenshtein(barcode, softclip_seq[i:i+bc_len]))
      if dist[-1] == 1: break  
    edit_dist = min(dist)
    bc_idx = dist.index(edit_dist)
    return [edit_dist, i_range[bc_idx]]
  
def search_pos_order(search_len, search_start, token_len):
  #Normally for reverse strand gene, search_len=55, start=23, token_len(~barcode)=16
  #  for plus strand gene, search_len=55, start=len(search_seq)-39=16
  min_pos=0; max_pos=search_len-token_len;
  search_range = [search_start]
  
  max_i = max(max_pos-token_len, min_pos+token_len)
  for i in range(1, max_i+1):
    if search_start+i <= max_pos:
      search_range.append(search_start+i)
    if search_start-i >= min_pos:
      search_range.append(search_start-i)
  return search_range

def find_r1_adapter(softclip_seq, adapter_bp, bc_pos, strand):
   if strand == 'plus':
     r1_seq_start = max(bc_pos-adapter_bp-2, 0)
     r1_seq_end = min(bc_pos+2, len(softclip_seq))
     pos_adjust = len(softclip_seq)  # Adjustment to return plus strand gene position as negative offset from end of softclip sequence
     adapter_seq = ILLUMINA_R1_FWD[-adapter_bp:]
     #print("plus:", bc_pos, adapter_seq, r1_seq_start, r1_seq_end)

   else:
     r1_seq_start = bc_pos+16
     r1_seq_end = min(r1_seq_start+adapter_bp+2, len(softclip_seq))
     pos_adjust = 0
     adapter_seq = ILLUMINA_R1_REV[0:adapter_bp]
     #print("minus:", bc_pos, adapter_seq, r1_seq_start, r1_seq_end)

   found_adapter_pos = [m.start(0) for m in re.finditer(adapter_seq, softclip_seq)]
   if len(found_adapter_pos) == 1:
     return [0, found_adapter_pos[0] - pos_adjust]
   elif len(found_adapter_pos) > 0:
     return [0, closest_pos_match(found_adapter_pos, r1_seq_start) - pos_adjust]

   else:
     if r1_seq_end > r1_seq_start + adapter_bp -1:
       max_i = max(r1_seq_end+1-adapter_bp, r1_seq_start) 
       #print("Searching for adapter sequence from:", r1_seq_start, "to:", r1_seq_end, "in softclip length:", len(softclip_seq), "BCpos:", bc_pos) 
       dist = [];  
       for i in range(r1_seq_start, max_i):
         dist.append(distance.levenshtein(adapter_seq, softclip_seq[i:i+adapter_bp]))
       edit_dist = min(dist)
       return [edit_dist, dist.index(edit_dist) + r1_seq_start - pos_adjust]
     
     else:
       return [-1, -1]

