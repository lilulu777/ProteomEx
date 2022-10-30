pacman::p_unload(pacman::p_loaded(), character.only = T)
library(tidyverse)
library(magrittr)

df_mat <- rio::import('ProteomEx_demo_protein_10%_clean.xlsx')
df_lbl <- rio::import('ProteomEx_protMatrix_demoAll_Label_20221028.csv')

df_merge <- df_mat %>%
  column_to_rownames('prot') %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('SampleID') %>%
  full_join(df_lbl, ., by = 'SampleID')
rio::export(df_merge, 'ProteomEx_demo_protein_10%_clean_addLabel_20221029.xlsx')
  




