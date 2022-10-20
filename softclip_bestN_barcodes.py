#!/usr/bin/env python
"""

:Author: Ji Research Group/Stanford Genome Technology Center
:Contact: sgrimes@stanford.edu
:Creation date: 03/24/2021
:Description: 

This script extracts soft clipped bases at beginning (FWD strand) or end (REV strand)
of read.  These sequences will subsequently be searched for expected single cell barcode sequences.
    
Revisions: 

- 03/26/2021	Import reusable methods from sc_barcodes
- 05/10/2021	Start barcode search from expected start position (vs beginning of sequence)
- 05/11/2021	Add nbest and threshold parameters
- 02/24/2022	Modify to include extra SC and aligned bases if UMI < 10 after barcode matching

"""
import argparse, sys, os, re, pysam, csv, gzip, string, distance
import numpy as np, pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from scipy.sparse import csr_matrix

import sc_barcodes as scb

script_name = os.path.basename(__file__)
print("Running ", script_name)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define internal modules                                                     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+# 
def parse_commandline():
  parser=argparse.ArgumentParser()
  #parser.add_argument('--bam', '-b', help='input bam file', type=argparse.FileType('r'), required=True)
  parser.add_argument('--bam', '-b', help='input bam file', type=str, required=True)
  parser.add_argument('--barcodes', '-c', help='cellranger barcodes file', type=str, required=True)
  parser.add_argument('--strand', '-s', help='gene strand', type=str, required=True, choices=['plus','minus','both'])
  parser.add_argument('--exonrds', '-x', help='reads with exon skipping pattern identified', type=str, required=False)
  parser.add_argument('--kmer_len', '-k', help='k-mer length', type=int, required=False, default=8)
  parser.add_argument('--offset', '-o', help='search offset from beg/end of soft clip', type=int, required=False, default=9)
  parser.add_argument('--nbest', '-n', help='max number of top cos-sim scores to evaluate', type=int, required=False, default=5)
  parser.add_argument('--threshold', '-t', help='cosine similarity min score threshold', type=int, required=False, default=0.1)
 
  args=parser.parse_args()
  print(args, file=sys.stderr)
  return args
  
def debug_samrd(samrd):
  strand = 'Rev' if samrd.is_reverse else 'Fwd'
  return [samrd.qname[-16:], strand, samrd.tid, samrd.pos]
  
def best_barcodes(string, barcodes, barcode_tfidf, vectorizer, nbest=5):
  best_barcodes = [['N',0]]
  barcode_seq_tfidf = vectorizer.transform([string])
  cos_sim = cosine_similarity(barcode_seq_tfidf, barcode_tfidf, dense_output=False)
  
  non_zero = [((i, j), cos_sim[i,j]) for i, j in zip(*cos_sim.nonzero())]
  nz_sorted = sorted(non_zero, key=lambda x: -x[1])
  idx_nbest = [x[0][1] for x in nz_sorted[0:nbest] if x[1] > MIN_SCORE]
  
  if len(idx_nbest) > 0:
    best_barcodes = zip([barcodes[i] for i in idx_nbest], [cos_sim[(0,i)] for i in idx_nbest])
  
  return best_barcodes

def format_bc_string(soft_clips, strand, pos, aligned_bases=None):
  if strand == 'fwd':  #positions will all be -ve offsets from end of sequence
    flank_umi = soft_clips[pos+16:min(pos+26, len(soft_clips))]
    if len(flank_umi) < 10 and aligned_bases:
      umi = ''.join([flank_umi, aligned_bases[0:10-len(flank_umi)]])
      #print("(fwd) orig umi:", flank_umi, "new umi: ", umi)
    else:
      umi = flank_umi

    r1_adapter = soft_clips[max(pos-22, 0):pos]
    barcode = soft_clips[pos:pos+16]
    return '|'.join([r1_adapter, barcode, umi])
  
  else:
    flank_umi = soft_clips[max(pos-10, 0):pos]
    if len(flank_umi) < 10 and aligned_bases:
      xstart = len(aligned_bases) - (10 - len(flank_umi))
      umi = ''.join([aligned_bases[xstart:], flank_umi])
      #print("(rev) orig umi:", flank_umi, "[", scb.reverse_complement(flank_umi), "]", "new umi: ", umi, "[", scb.reverse_complement(umi), "]")
    else:
      umi = flank_umi

    barcode = soft_clips[pos:pos+16]
    r1_adapter = soft_clips[pos+16:min(pos+22, len(soft_clips))]
    return '|'.join([scb.reverse_complement(r1_adapter), scb.reverse_complement(barcode), scb.reverse_complement(umi)])

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
args = parse_commandline()

#Use SEARCH_OFFSET to ignore part of 13bp 10X adapter, and 10bp UMI (9 gives generous allowance for indels and partial alignment of 10x adapter)
SEARCH_OFFSET = args.offset
MAX_SEARCH_LEN = 55-SEARCH_OFFSET
KMER_LEN = args.kmer_len
NBEST = args.nbest
MIN_SCORE = args.threshold
ADAPTER_BP = 7

sam_input = pysam.Samfile(args.bam,'rb') if args.bam[-3:] == 'bam' else pysam.Samfile(args.bam,'r')
sam_fname = os.path.basename(args.bam)

if args.barcodes.endswith('gz'):
  barcode_input = gzip.open(args.barcodes, 'r')
  barcodes = [line[:16].decode('ascii') for line in barcode_input]
else:
  barcode_input = open(args.barcodes, 'r')
  barcodes = [line[:16] for line in barcode_input]

xskip_fn = args.exonrds

try:
  out_fn = sam_fname[0:-4] + '.softclip.bestN.txt'
  out_file = open(out_fn, 'w')
  out_csv = csv.writer(out_file, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
except:
  print("Unable to open text file for output: ", out_fn)
  sys.exit(1)
  
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Vectorize barcodes                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#https://bergvca.github.io/2017/10/14/super-fast-string-matching.html
def ngrams(string, kmer_len=KMER_LEN):
  ngrams = zip(*[string[i:] for i in range(kmer_len)])
  return [''.join(ngram) for ngram in ngrams]
  
#build kmer dictionary of all barcode seqs (in both forward and reverse orientation)
vectorizer_fwd = CountVectorizer(min_df=1, analyzer=ngrams)
fwd_barcode_tfidf = vectorizer_fwd.fit_transform(barcodes)

vectorizer_rev = CountVectorizer(min_df=1, analyzer=ngrams)
rev_barcode_tfidf = vectorizer_rev.fit_transform([scb.reverse_complement(barcode) for barcode in barcodes])

if args.exonrds: 
  xskip_rdnames = pd.read_csv(xskip_fn, sep='\t', header=None, index_col=0)
   
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Read sam file and check for soft-clips at beginning or end of read          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Expected sequence (FWD), soft clips indicated by []
#   [Illumina R1 adapter][cell barcode][UMI][10X internal adapter] cDNA [polyA][adapter]
#
# Expected sequence (REV), soft clips indicated by []
#   [adapter][polyT] cDNA [10X internal adapter][UMI][barcode][Illumina R1 adapter]
i = 0; tot_rds = 0; adapter_flank_bp = 6;

out_csv.writerow(['rd_name', 'exon_skip', 'strand', 'barcode', 'score', 'dist', 'pos', 'pos_offset', 'adapter|BC|UMI', 
                  'r1_dist', 'r1_pos', 'search_len', 'sc_orientation'])
sc_3or5 = '5prime' if args.strand == 'plus' else '3prime'

for samrd in sam_input.fetch(until_eof=True):
  i += 1
  
  if samrd.is_secondary:
    continue

  if args.exonrds and samrd.qname in xskip_rdnames.index:
    xskip_pattern = xskip_rdnames.loc[samrd.qname,1]
  elif args.exonrds:
    continue
  else:
    xskip_pattern = None

  align_strand = 'minus' if samrd.is_reverse else 'plus'

  tot_rds += 1
  sam_flds = debug_samrd(samrd)
  if i < 5:
    print("Checking for soft clips: {0}".format(sam_flds))
  
  soft_clips = scb.extract_softclips(samrd)
  aligned_10 = scb.first_aligned_bases(samrd, n=10)
 
  barcode_scores_fwd = [['N', 0]]; barcode_scores_rev = [['N', 0]];
  search_seq_fwd = ''; search_seq_rev = '';

  if args.strand in ['both', 'plus']:
    sc_5prime_len = len(soft_clips['fwd'])
    if sc_5prime_len > 16+SEARCH_OFFSET+1:    
      #Working backwards from end of soft clipped sequence (s/b 10X adapter, then UMI, then cell barcode) - determine start position for cell barcode search
      i_start = max(sc_5prime_len-MAX_SEARCH_LEN-SEARCH_OFFSET, 0)
      search_seq_fwd = soft_clips['fwd'][i_start:-SEARCH_OFFSET] if SEARCH_OFFSET > 0 else soft_clips['fwd'][i_start:]
      extra_sc_fwd = soft_clips['fwd'][-SEARCH_OFFSET:]
      barcode_scores_fwd = best_barcodes(search_seq_fwd, barcodes, fwd_barcode_tfidf, vectorizer_fwd, NBEST)
          
  if args.strand in ['both', 'minus']: 
    sc_3prime_len = len(soft_clips['rev'])
    if sc_3prime_len > 16+SEARCH_OFFSET+1:     
      i_end = min(MAX_SEARCH_LEN+SEARCH_OFFSET, sc_3prime_len)  
      search_seq_rev = soft_clips['rev'][SEARCH_OFFSET:i_end]
      extra_sc_rev = soft_clips['rev'][0:SEARCH_OFFSET]
      barcode_scores_rev = best_barcodes(search_seq_rev, barcodes, rev_barcode_tfidf, vectorizer_rev, NBEST)
  
  for bc_score in barcode_scores_fwd:
    if bc_score[0] != "N":
      search_seq = search_seq_fwd
      search_start = max(0, len(search_seq)-39)  #TODO: Should adjust where edit distance searching starts, based on SEARCH_OFFSET
      [dist, pos] = scb.calc_edit_distance(bc_score[0], search_seq, search_start)
      [r1_dist, r1_pos] = scb.find_r1_adapter(search_seq, ADAPTER_BP, pos, args.strand)
      barcode_with_flanking = format_bc_string(search_seq, 'fwd', pos, "".join([extra_sc_fwd, aligned_10['fwd']]))
      bc_pos = pos-len(search_seq)
      bc_offset = abs(len(search_seq)-39 - pos)
      out_csv.writerow([samrd.qname, xskip_pattern, align_strand, bc_score[0], bc_score[1], dist, bc_pos, bc_offset, 
                        barcode_with_flanking, r1_dist, r1_pos, len(search_seq), 'left'])
      if dist <= 1: break

  for bc_score in barcode_scores_rev:
    if bc_score[0] != "N":
      search_seq = search_seq_rev
      search_start = min(len(search_seq)-16, 23) #TODO: Should adjust where edit distance searching starts, based on SEARCH_OFFSET
      [dist, pos] = scb.calc_edit_distance(scb.reverse_complement(bc_score[0]), search_seq, search_start)
      [r1_dist, r1_pos] = scb.find_r1_adapter(search_seq, ADAPTER_BP, pos, args.strand)
      barcode_with_flanking = format_bc_string(search_seq, 'rev', pos, "".join([aligned_10['rev'],extra_sc_rev]))
      bc_pos = pos
      bc_offset = abs(23 - pos)  
      out_csv.writerow([samrd.qname, xskip_pattern, align_strand, bc_score[0], bc_score[1], dist, bc_pos, bc_offset, 
                        barcode_with_flanking, r1_dist, r1_pos, len(search_seq), 'right'])
      if dist <= 1: break

print(i, "sam records read")
print("Evaluated", tot_rds, "primary (full transcript) alignments")

for fh in [sam_input, barcode_input, out_file]:
  fh.close()
