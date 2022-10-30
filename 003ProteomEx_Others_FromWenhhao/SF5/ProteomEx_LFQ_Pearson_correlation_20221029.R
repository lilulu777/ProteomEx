require(magrittr)
require(dplyr)
require(stringr)
library(tidyverse)

rm(list = ls())
# setwd('//172.16.13.114/share/members/jiangwenhao/lilu/ProteomEx/LFQ_correlation')

rlt_vec <- list.files('//172.16.13.136/ProteomEx/ProteomEx_LFQ_20220606/para_search/ProteomEx_LFQ', pattern = '^protein\\.tsv$', full.names = T, recursive = T)
filename_vec <- stringr::str_extract(rlt_vec, 'ProteomEx_[A-Za-z0-9_]+?_Slot')

df_ls <- lapply(rlt_vec, function(e){
  df_tmp <- data.table::fread(e, sep = '\t', data.table = F)
  df_ret <- df_tmp %>%
    select(Protein, `Unique Intensity`) %>%
    filter(!grepl('^CON', Protein))
  colnames(df_ret)[2] %<>% paste0(., '_', stringr::str_extract(e, 'ProteomEx_[A-Za-z0-9_]+?_Slot'))
  return(df_ret)
})

df_merge <- Reduce(function(x, y){
  merge(x, y, by = 'Protein', all = T)
}, df_ls)
df_merge[df_merge == 0] <- NA
pm <- df_merge[, -1]
rownames(pm) <- df_merge$Protein
colnames(pm) %<>% stringr::str_replace('Unique Intensity_ProteomEx_', '') %>% stringr::str_replace('_Slot', '')

# 列排序 ---------------------------------------------------------------------
pm %<>% select(c1, c2, c4, c3, b1:b4, G1_1_dz, G1_3_dz, G1_4_dz, G1_1_lilu, G1_2_lilu, G1_3_lilu, G1_4_lilu, a1:a4, everything())


# correlation -------------------------------------------------------------
require(corrplot)
cors <- cor(log2(pm), use = "pairwise.complete.obs", method = 'pearson')
uptri <- cors[upper.tri(cors)]
median(uptri) # 0.9650217
min(uptri) # 0.6258299
#hist(uptri)
pdf('ProteomEx_LFQ_Pearson_correlation_23files_myrank_20220616.pdf', width = 22, height = 22)
# corrplot(cors, order ="AOE", type = "upper", diag = F, addCoef.col = "grey", cl.cex = 2) # 太多了看不清
corrplot(cors, method = 'square', cl.cex = 2, tl.pos = 'lt', tl.col = 'black', tl.cex = 2.2, addCoef.col = 'grey', number.cex = 2)
dev.off()



# CV ----------------------------------------------------------------------
groups <- colnames(pm) %>% str_sub(1, 1) %>% unique
df_cv <- lapply(groups, function(e){
  mat_tmp <- pm[, grep(e, colnames(pm))]
  cv_tmp <- apply(mat_tmp, 1, function(x){sd(x, na.rm = T) / mean(x, na.rm = T)})
  return(cv_tmp)
}) %>% do.call(cbind, .)
colnames(df_cv) <- groups
df_cv <- cbind(data.frame(prot = rownames(df_cv)), df_cv)
df_vio <- reshape2::melt(df_cv, 'prot', variable = 'Group', value.name = 'CV')
df_vio %<>% tidyr::drop_na('CV')
df_fiv <- apply(df_cv[, -1], 2, fivenum, na.rm = T)
rownames(df_fiv) <- c('min', '1st quartile', 'median', '3rd quartile', 'max')
df_fiv <- apply(df_fiv, c(1, 2), function(s) as.numeric(sprintf('%0.2f', s)))
df_fiv

# Number n ----------------------------------------------------------------------
ht_group <- c(a = 'ProteomEx',
              b = 'PCT',
              c = 'In-solution',
              G = 'proExM-MS',
              S = 'ProteomEx-2'
              )

protn <- df_cv[, -1] %>%
  apply(2, function(col){
    c(protein_number = sum(!is.na(col)))
  })
names(protn) %<>% ht_group[.]
protn %>% t() %>% t() %>% as.data.frame() %>% setNames('protein_number') %>% rownames_to_column('Group') %>% rio::export('ProteomEx_SF6_correlation_protein_number_20221029.xlsx')
  
library(ggplot2)
med <- median(df_vio$CV)
p <- ggplot(df_vio, aes(x = Group, y = CV, fill = Group))+
  geom_violin() +
  theme(legend.position = 'none',
        #legend.text = element_text(size = 20,color = 'black'),legend.position = 'right',
        #legend.title = element_text(size = 20,color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  theme(panel.grid = element_blank())+
  labs(
    #title = '',
    #subtitle = '',
    #caption = '',
    #tag = '',
    #x = '',
    y = 'Protein intensity CV')+
  theme(axis.text = element_text(size = 20,color = 'black'))+
  theme(axis.text.x = element_text(size = 20,color = 'black'))+
  theme(plot.subtitle = element_text(size = 25, hjust = 0, color = 'black'))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(size = 22, hjust = 0.5, color = 'black'))+
  #scale_y_continuous(name = 'log2Ratio', breaks = seq(-3, 3, 1)) +
  scale_fill_brewer(palette = 'Pastel2')+
  # scale_fill_manual(values = c('#E05B28'))
  stat_summary(fun = mean, geom = 'point', size = 2, color = 'red')

p <- p + geom_boxplot(width = 0.1) + stat_summary(fun = mean, geom = 'point', size = 2, color = 'red')
#+annotate("text", x = 1, y = max(df_vio$CV, na.rm = T) * 1.05, parse = T, label = stringr::str_c('CV-Median: ', sprintf('%.3f', med)), size = 6, color = 'darkred')
ggsave('20220616ProteomEx_CV_raw.pdf', p, width = 6.2, height = 3.64)



# CV after log2 ----------------------------------------------------------
pm <- log2(pm)
groups <- colnames(pm) %>% str_sub(1, 1) %>% unique
df_cv <- lapply(groups, function(e){
  mat_tmp <- pm[, grep(e, colnames(pm))]
  cv_tmp <- apply(mat_tmp, 1, function(x){sd(x, na.rm = T) / mean(x, na.rm = T)})
  return(cv_tmp)
}) %>% do.call(cbind, .)
colnames(df_cv) <- groups
df_cv <- cbind(data.frame(prot = rownames(df_cv)), df_cv)
df_vio <- reshape2::melt(df_cv, 'prot', variable = 'Group', value.name = 'CV')
df_vio %<>% tidyr::drop_na('CV')

library(ggplot2)
med <- median(df_vio$CV)

df_med <- df_vio %>%
  group_by(Group) %>%
  summarise(median = median(CV))


p <- ggplot(df_vio, aes(x = Group, y = CV, fill = Group))+
  geom_violin() +
  theme(legend.position = 'none',
        #legend.text = element_text(size = 20,color = 'black'),legend.position = 'right',
        #legend.title = element_text(size = 20,color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  theme(panel.grid = element_blank())+
  labs(
    #title = '',
    #subtitle = '',
    #caption = '',
    #tag = '',
    #x = '',
    y = 'Protein intensity CV')+
  theme(axis.text = element_text(size = 20,color = 'black'))+
  theme(axis.text.x = element_text(size = 20,color = 'black'))+
  theme(plot.subtitle = element_text(size = 25, hjust = 0, color = 'black'))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(size = 22, hjust = 0.5, color = 'black'))+
  #scale_y_continuous(name = 'log2Ratio', breaks = seq(-3, 3, 1)) +
  scale_fill_brewer(palette = 'Pastel2')+
  # scale_fill_manual(values = c('#E05B28'))
  stat_summary(fun = mean, geom = 'point', size = 2, color = 'red')

p <- p + geom_boxplot(width = 0.1) + stat_summary(fun = mean, geom = 'point', size = 2, color = 'red')
  #+annotate("text", x = 1, y = max(df_vio$CV, na.rm = T) * 1.05, parse = T, label = stringr::str_c('CV-Median: ', sprintf('%.3f', med)), size = 6, color = 'darkred')
ggsave('20220616ProteomEx_CV_log2.pdf', p, width = 6.2, height = 3.64)

