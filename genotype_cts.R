library(dplyr)
library(tidyr)
library(data.table)
library(stringdist)
firstInList <- function(x) { sapply(x, "[[", 1)[1]}

# Guide#s for for loop
targets <- read.csv('NGG_sgRNA_coordination.txt', sep='\t')
sample_prefix <- 'Sample_Name'
guides <- targets$Target_site

#Barcodes by ONT readname
bc_info_file <- paste0(sample_prefix, '.barcode_match.tsv')
bc_info <- fread(bc_info_file, sep='\t', header=T, stringsAsFactors=F)

rd_barcodes <- bc_info %>% filter(dist==0) %>%
  select(rd_name, barcode, dist, "adapter|BC|UMI") %>%
  separate("adapter|BC|UMI", c('r1_flank', 'barcode_seq', 'umi'))

Ref_seq <- targets$Ref_seq

for (sgrna in guides) {
  print('start analysis of guideRNA')
  print(sgrna)
  #Merge ONT readname/genotypes with cell barcode/UMIs
  in_fn <- paste0(sample_prefix,'.TP53_sgrna', sgrna, '.be_genotypes.txt')
  tp53GT <- fread(in_fn, header=F, sep='\t', stringsAsFactors=F)
  colnames(tp53GT) <- c('rd_name', 'genotype')
  
  #Remove unexpected bases
  ref_seq <- Ref_seq[sgrna]
  target_base <- substr(targets[sgrna,]$Targeted.Base,1,1)
  grep_pattern = gsub(target_base, "[ACGT]", ref_seq, fixed=T)
  tp53GT$gt_match <- apply(tp53GT, 1, FUN=function(x) grep(grep_pattern, x['genotype']))
  tp53GT <- tp53GT %>% filter(gt_match == 1) %>% select(-gt_match)
  
  bc_genotype <- inner_join(x=tp53GT, y=rd_barcodes)
  
  #UMI de-duplication and cell barcode genotyping
  bc_umi <- bc_genotype %>% 
    select(barcode, r1_flank, barcode_seq, umi, genotype, dist)
  out_fn <- paste0(sample_prefix,'.TP53_sgrna', sgrna, '.bc_umi.tsv')
  write.table(bc_umi, file=out_fn, sep='\t', row.names=F, quote=F)
    
  print('Starting bc_dedup')
  
  bc_temp <- bc_umi %>% group_by(barcode, umi) %>% count(genotype) %>% filter(n>9) %>% arrange(desc(n)) %>% 
    top_n(2) %>% filter(row_number() < 3) %>%
    mutate(wt_mut = ifelse(genotype == ref_seq, 'WT', 'Mut'), pct=prop.table(n))
  
  bc_dedup <- bc_temp %>%
    pivot_wider(names_from = wt_mut, values_from=c('genotype', 'n', 'pct'), values_fn=firstInList) %>%  
    mutate(final_GT = ifelse((is.na(pct_WT) | pct_WT <= 0.85), genotype_Mut, genotype_WT)) %>%
    select(barcode, umi, final_GT) %>% rename(genotype = final_GT)
  
  print('Starting bc_consolidation')
  
  barcode_list <- unique(bc_dedup$barcode)
  bc_umi_consolidate <- data.frame(matrix(nrow = 0, ncol = length(colnames(bc_dedup))))
  colnames(bc_umi_consolidate) <- colnames(bc_dedup)
  for (bc in barcode_list){
    umi_idx = 1
    tmp_df <- bc_dedup %>% filter(barcode==bc)
    while (umi_idx < nrow(tmp_df)) {
      tmp_df <- tmp_df %>% 
        mutate(row_id = row_number(), ED=stringdist(tmp_df$umi[umi_idx], umi, method='lv')) %>%
        filter(row_id <= umi_idx, ED==0 | ED>2) %>%
        select(-c(ED, row_id))
      umi_idx = umi_idx + 1
    }
    bc_umi_consolidate <- rbind(bc_umi_consolidate, tmp_df)
  }

  genotype_cts <- bc_umi_consolidate %>% group_by(barcode, genotype) %>% 
                               summarize(nr_rds = n()) %>% arrange(desc(nr_rds), .by_group=T)
  out_fn1 <- paste0(sample_prefix,'.TP53_sgrna', sgrna, '.genotype_cts.tsv')
  write.table(genotype_cts, file=out_fn1, sep='\t', row.names=F, quote=F)
  
  #Summary genotype counts for sgRNA
  out_fn2 <- paste0('sgrna', sgrna, '.genotype_summary.tsv')
  summary_cts <- genotype_cts %>% group_by(genotype) %>% summarize(tot_rds = sum(nr_rds)) %>% arrange(desc(tot_rds))
  write.table(summary_cts, file=out_fn2, sep='\t', row.names=F, quote=F)
}


