library(tidyverse)
library(magrittr)
# require(ggsignif)

rm(list = ls())

# my_colors <- c(`In-solution` = '#064F89', 
#                PCT = '#E6B429',
#                `proExM-MS` = '#8FC31F', 
#                ProteomEx = '#DD704A',
#                `ProteomEx-2` = '#DD704A')
# 
# group_rename <- c('In-solution', 'PCT', 'proExM-MS', 'ProteomEx')
# names(group_rename) <- c('c', 'b', 'g', 'a')


files_AGE <- list.files('//172.16.13.136/ProteomEx/ProteomEx_LFQ_20220606/20220806_AGE', pattern = '^peptide\\.tsv$', full.names = T, recursive = T)
files_NSA <- list.files('//172.16.13.136/ProteomEx/ProteomEx_LFQ_20220606/20220806_NSA', pattern = '^peptide\\.tsv$', full.names = T, recursive = T)

rlt_ls_AGE <- lapply(files_AGE, function(f){
  df <- read.delim(f, check.names = F)
  intrested_sites <- c('K', 'R', 'N', 'Q')
  data.frame(
    Ratio = sapply(intrested_sites, function(pattern){
      sum(str_detect(df$Peptide, pattern))
    }) / nrow(df),
    `Peptide number with this aa` = sapply(intrested_sites, function(pattern){
      sum(str_detect(df$Peptide, pattern))
    })
  )
})

rlt_ls_NSA <- lapply(files_NSA, function(f){
  df <- read.delim(f, check.names = F)
  intrested_sites <- c('K', 'R', 'N', 'Q')
  data.frame(
    Ratio = sapply(intrested_sites, function(pattern){
      sum(str_detect(df$Peptide, pattern))
    }) / nrow(df),
    `Peptide number with this aa` = sapply(intrested_sites, function(pattern){
      sum(str_detect(df$Peptide, pattern))
    })
  )
})

df_rlt_AGE <- Reduce(cbind, rlt_ls_AGE)
df_rlt_NSA <- Reduce(cbind, rlt_ls_NSA)

colnames(df_rlt_AGE) <- str_split(files_AGE, '/') %>% lapply(tail, 2) %>% sapply(head, 1) %>% rep(each = 2) %>% str_c(c('ratio', 'number'), sep = '_')
colnames(df_rlt_NSA) <- str_split(files_NSA, '/') %>% lapply(tail, 2) %>% sapply(head, 1) %>% rep(each = 2) %>% str_c(c('ratio', 'number'), sep = '_')

df_rlt <- cbind(df_rlt_AGE, df_rlt_NSA) %>% rownames_to_column('Amino acid')
rio::export(df_rlt, 'ProteomEx_AGE-NSA_KRQN_ratio_SF11B_20220811.xlsx')

