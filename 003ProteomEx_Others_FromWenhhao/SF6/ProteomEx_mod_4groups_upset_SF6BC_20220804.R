require(magrittr)
library(tidyverse)

df <- read.delim('//172.16.13.136/ProteomEx/ProteomEx_LFQ_20220606/20220804ProteomEx_PTM_canonical_variable_modification/global.modsummary.tsv', sep = '\t', check.names = F, stringsAsFactors = F)
df1 <- df %>% 
  select(Modification, a_percent_PSMs:g_percent_PSMs) %>%
  dplyr::rename('ProteomEx' = a_percent_PSMs,
         'PCT' = b_percent_PSMs,
         'In-solution' = c_percent_PSMs,
         'proExM-MS' = g_percent_PSMs) %>%
  pivot_longer(-Modification, names_to = 'group', values_to = 'percent')
  # reshape2::melt(id = 'Modification', var = 'group', value.name = 'percent')

# library(plyr)
# df2 <- plyr::ddply(df1, 'group', function(dfsub){
#   dfsub %>% arrange(desc(percent)) %>%
#     slice(1:50)
# })
# 
# df2 %>% dplyr::count(Modification)

# upset
library(UpSetR)
# listInput <- list(
#   one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), 
#   two = c(1, 2, 4, 5, 10), 
#   three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))
df3 <- df1 %>% reshape2::dcast(Modification ~ group) %>%
  select(Modification, 4, 3, 2, 5) # 长转宽
ls3_upset <- apply(df3[, -1], 2, function(col){
  df3$Modification[col != 0]
}, simplify = F)

pdf('ProteomEx_mod_4groups_upset_SF6B_20220805.pdf', width = 8, height = 5)
print(UpSetR::upset(UpSetR::fromList(ls3_upset), order.by = 'freq'))
dev.off()


# heatmap -----------------------------------------------------------------
require(RColorBrewer)
mat <- df3 %>% filter(Modification != 'None') %>% column_to_rownames('Modification')
pheatmap::pheatmap(mat, cluster_rows = T,
                   color = c(brewer.pal(11,"RdYlBu")[9:7],brewer.pal(11,"RdYlBu")[4:2]))

mat_log <- mat * 1000
sort(unique(unlist(mat_log))) #  0    1    2    3    4    5    6    7   .....
mat_log[mat_log == 0] <- 1
mat_log %<>% log10
pheatmap::pheatmap(mat_log, color = c(brewer.pal(11,"RdYlBu")[9:7],brewer.pal(11,"RdYlBu")[4:2]), fontsize_col = 8,
                   # annotation_col = ann_col, annotation_row = ann_row,
                   # gaps_col = cumsum(table(ann_col$tissue))[-length(cumsum(table(ann_col$tissue)))],
                   # gaps_row = cumsum(table(row.tissue))[-length(cumsum(table(row.tissue)))],
                   # fontsize_row=3,
                   # annotation_colors = ann_colors,
                   na_col = "grey60",#scale = "row",
                   cluster_rows = T, cluster_cols = T,show_rownames = F, show_colnames = T, 
                   filename = "ProteomEx_mod_4groups_heatmap_SF6C_20220805.pdf",width=8,height=8)

get_label <- function(x){(10 ^ x) / 1000}
get_label(0:3) # 0.001 0.010 0.100 1.000 (% ; in fact)




# check interesting points
only_pct_insol <- intersect(ls3_upset$PCT, ls3_upset$`In-solution`) %>%
  setdiff(ls3_upset$`proExM-MS`) %>%
  setdiff(ls3_upset$ProteomEx)

df %>% filter(Modification %in% only_pct_insol)






# 用不上了 --------------------------------------------------------------------


# 
# mat.scale <- apply(mat, 1, scale) %>% t() %>% data.frame()
# mat.scale[is.na(mat.scale)] <- 0
# colnames(mat.scale) <- colnames(df3)[-1]
# pheatmap::pheatmap(mat.scale, cluster_rows = T,
#                    color = c(brewer.pal(11,"RdYlBu")[9:7],brewer.pal(11,"RdYlBu")[4:2]))



# 
# pheatmap(Protein_Matrix_DIANN_order,
#          breaks = bk, border_color = F,
#          scale = "row",
#          gaps_col = cumsum(table(anno_col$Group))[-length(cumsum(table(anno_col$Group)))],
#          gaps_row = cumsum(table(anno_row$Cluster))[-length(cumsum(table(anno_row$Cluster)))],
#          fontsize = 9,
#          color = my_colors,
#          annotation_col = anno_col,
#          annotation_row = anno_row,
#          annotation_colors = ann_colors,
#          show_rownames = F, show_colnames = T, 
#          width = 8,height = 8,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          filename = 'MCAO_heatmap_gaps_CM012_3clusters_20220715.pdf'
# )




# proteomex_prot <- ls3_upset$ProteomEx
# insolution_pct_prot <- intersect(ls3_upset$`In-solution`, ls3_upset$PCT)
# proteomex_unique_prot <- setdiff(proteomex_prot,
#                                  Reduce(union,
#                                         list(
#                                           ls3_upset$`In-solution`,
#                                           ls3_upset$PCT,
#                                           ls3_upset$`proExM-MS`
#                                         ))
#                                  )
# 
# insolution_pct_group_unique_prot <- setdiff(insolution_pct_prot,
#                                             union(ls3_upset$ProteomEx,
#                                                   ls3_upset$`proExM-MS`))




# # rank
# df_percent <- df[-1, ] %>% # remove None modification
#   select(a_percent_PSMs:g_percent_PSMs) %>%
#   dplyr::rename('ProteomEx' = a_percent_PSMs,
#                 'PCT' = b_percent_PSMs,
#                 'In-solution' = c_percent_PSMs,
#                 'proExM-MS' = g_percent_PSMs) %>%
#   apply(2, sort, decreasing = T) %>%
#   as.data.frame %>%
#   tibble::rownames_to_column('rank') %>%
#   reshape2::melt(id = 'rank', var = 'group', value.name = 'percent')
# df_percent$rank %<>% as.numeric()
# 
# library(ggplot2)
# p <- ggplot(df_percent) +
#   aes(x = rank, y = percent, fill = group, colour = group) +
#   geom_point(shape = "circle", size = 1.5) +
#   scale_fill_hue(direction = 1) +
#   scale_color_hue(direction = 1) +
#   theme_classic()+
#   theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+
#   theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"))+
#   theme(plot.title = element_text(size = 20, face = "bold"))+
#   labs(x="Rank",y="Percent (%)")+
#   theme(legend.position='right')
# 
# 
# # 按 ProteomEx 顺序
# df_percent2 <- df[-1, ] %>% # remove None modification
#   select(a_percent_PSMs:g_percent_PSMs) %>%
#   dplyr::rename('ProteomEx' = a_percent_PSMs,
#                 'PCT' = b_percent_PSMs,
#                 'In-solution' = c_percent_PSMs,
#                 'proExM-MS' = g_percent_PSMs) %>%
#   arrange(desc(ProteomEx)) %>%
#   tibble::rownames_to_column('rank') %>%
#   reshape2::melt(id = 'rank', var = 'group', value.name = 'percent')
# df_percent2$rank %<>% as.numeric()
# 
# p2 <- ggplot(df_percent2) +
#   aes(x = rank, y = percent, fill = group, colour = group) +
#   geom_point(shape = "circle", size = 1.5) +
#   scale_fill_hue(direction = 1) +
#   scale_color_hue(direction = 1) +
#   theme_classic()+
#   theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+
#   theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"))+
#   theme(plot.title = element_text(size = 20, face = "bold"))+
#   labs(x="Modification",y="Percent (%)")+
#   theme(legend.position='right')
# 
# pdf('modification_4groups_top50_rank.pdf', width = 12, height = 7)
# print(p2)
# dev.off()


# # check interesting points ------------------------------------------------
# 
# require(tidyr)
# # wider table ---------------
# tb1 <- df_percent2 %>%
#   pivot_wider(names_from = group, values_from = percent) %>% arrange(desc(ProteomEx))
# tb1[is.na(tb1)] <- 0
# 
# # hash table ---------------
# ht <- df1 %>% filter(group == 'ProteomEx') %>%
#   arrange(desc(percent))
# ht_rank_prot <- ht$Modification
# ht_rank_prot <- ht_rank_prot[-1]
# names(ht_rank_prot) <- seq_along(ht_rank_prot)
# 
# # for preExM-MS ---------------
# tb1 %>%
#   filter(rank > 20, `proExM-MS` > 0.5, `In-solution` < 0.5)
# 
# p2 + geom_hline(yintercept=0.5, 
#                 color = "red", size=0.2)+
#   geom_vline(xintercept=c(35, 37, 47, 75, 93), 
#              color = "red", size=0.2)
# 
# ht_rank_prot[c(35, 37, 47, 75, 93)]
# # df3 %>% arrange(desc(ProteomEx)) %>% slice(-1) %>% slice(c(35, 37, 47, 75, 93)) %>% pull(Modification)
# 
# 
# # 35 
# # "Unannotated mass-shift 2.9744" 
# # 37 
# # "Iodoacetamide derivative/Addition of Glycine/Addition of G" 
# # 47 
# # "Homoserine lactone/Prompt loss of side chain from oxidised Met" 
# # 75 
# # "amidination of lysines or N-terminal amines with methyl acetimidate" 
# # 93 
# # "Ammonium-quenched monolink of BS2-G crosslinker" 
# 
# 
# 
# 
# # for PCT ---------------
# tb1 %>%
#   filter(rank > 50, PCT > 0.5)
# 
# p2 + geom_hline(yintercept=0.5, 
#                 color = "red", size=0.2)+
#   geom_vline(xintercept=c(72, 78), 
#              color = "red", size=0.2)
# 
# ht_rank_prot[c(72, 78)]
# 
# 
# 
# # 72                                78 
# # "Unannotated mass-shift -32.0250" "Unannotated mass-shift -59.0376" 
# 



