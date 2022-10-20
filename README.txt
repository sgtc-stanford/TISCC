# Running of scripts

1) Determine transcript genotypes per long read, for a region provided on command line
   Output: <output_prefix>.be_genotypes.txt
      
   python ${CODE_DIR}/transcript_genotypes.py <sample>.bam chr1:100000-100050 <output_prefix>   

2) Determine up to top5 barcode matches per long read, using vectorization and cosine similarity metric
   Required input files are bam file, and barcode 'whitelist' from short-reads.
   If transcript isoform file is specified, also merge in the transcript isoform per read, and exclude
     any reads which do not span the exons of interest.
   Output: <sample>.softclip.bestN.txt

   python softclip_bestN_barcodes.py -b <sample>.bam -c <sample>.barcodes.tsv.gz -s plus
   
   Parameters available with --help, or see below.
   '--bam', '-b', help='input bam file', type=str, required=True
   '--barcodes', '-c', help='cellranger barcodes file', type=str, required=True
   '--strand', '-s', help='gene strand', type=str, required=True, choices=['plus','minus','both']
   '--exonrds', '-x', help='reads with exon skipping pattern identified', type=str, required=False
   '--kmer_len', '-k', help='k-mer length', type=int, required=False, default=8
   '--nbest', '-n', help='number of best matches to evaluate', type=int, required=False, default=5

3) Pick best of top5 barcode matches, using edit distance as primary criterion (discard reads with edit distance > 4)
   Output: <sample>.barcode_match.tsv

   python select_best_barcode.py -i <sample>.softclip.bestN.txt

4) Final filtering of barcode matches.  Keep all barcode matches with edit distance < 3, and keep edit distance 3,4
     only if cosine similarity metric > 0.15.  Reorder columns and fix column heading for later processing.
   Output: <sample>.barcode_info.txt

   awk -F'\t' '{if(NR == 1 || ($6 < 3 || $5 > 0.15)) print $4, $2, $5, $3, $6, $7, $8, $9}' \
     <sample>.barcode_match.tsv > <sample>.barcode_info.txt
   sed -i '1 s/^.*$/best_barcode exon_skip max_score strand edit_dist start_pos barcode+flanking search_len/' <sample>.barcode_info.txt

5) Finalize genotyping for single-cell.  Cell-barcode per readID from '.barcode_match.tsv' file and genotype per readID 
     from '.be_genotypes.txt' are merged, then data is cleaned, deduplicated and counted.
   Input:  NGG_sgRNA_coordinates.txt or NG_sgRNA_coordinates.txt 
   Output: <sgRNA>.genotype_cts.tsv - single-cell genotype results for each sgRNA

   Rscript: genotype_cts.R

# License
Software in this repository is distributed according to the terms of the MIT license, as provided in the LICENSE file.


   
