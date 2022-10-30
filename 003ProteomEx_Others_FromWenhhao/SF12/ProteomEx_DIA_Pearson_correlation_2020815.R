require(magrittr)
require(tidyverse)


rm(list = ls())
df <- rio::import('20220810_proteins_LR_P2P.txt') %>% column_to_rownames('prot') # log2 transformed
# df %>% summarise_at(vars(everything()), function(x) {sum(!is.na(x))}) %>% t()
# df %>% colnames()
# df %>% select(matches('(Gel\\d_[^1]mm)|(PCT\\d\\.\\d_\\d)')) %>% colnames()

# remove Gel_1mm
df1 <- df %>% select(matches('(Gel\\d_[^1]mm)|(PCT\\d\\.\\d_\\d)')) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Filename')

df1_labeled <- df1 %>%
  add_column(
    Group = str_extract(df1$Filename, '(Gel)|(PCT)'),
    Amount = str_extract(df1$Filename, '(\\dmm)|(\\d\\.\\d)'),
    .before = 2
  ) %>%
  arrange(Group, Amount, Filename)

df_label <- df1_labeled %>% select(1:3)
pm <- df1_labeled %>%
  select(-(2:3)) %>%
  column_to_rownames('Filename') %>%
  t() %>%
  as.data.frame()



# correlation -------------------------------------------------------------
# require(corrplot)
cors <- cor(pm, use = "pairwise.complete.obs", method = 'pearson')
uptri <- cors[upper.tri(cors)]
median(uptri) # 0.7349846
min(uptri) # 0.5098424
#hist(uptri)
pdf('ProteomEx_DIA_Pearson_correlation_21files_20220815.pdf', width = 22, height = 22)
# corrplot(cors, order ="AOE", type = "upper", diag = F, addCoef.col = "grey", cl.cex = 2) # 太多了看不清
corrplot::corrplot(cors, method = 'square', cl.cex = 2, tl.pos = 'lt', tl.col = 'black', tl.cex = 2.2, addCoef.col = 'grey', number.cex = 2)
dev.off()

#Gel only
pm_gel <- df1_labeled %>%
  filter(str_detect(Group, 'Gel')) %>%
  select(-(2:3)) %>%
  column_to_rownames('Filename') %>%
  t() %>%
  as.data.frame()
cors <- cor(pm_gel, use = "pairwise.complete.obs", method = 'pearson')
uptri <- cors[upper.tri(cors)]
median(uptri) # 0.9309235
min(uptri) # 0.8859562
#hist(uptri)
pdf('ProteomEx_DIA_Pearson_correlation_Gel_12files_20220815.pdf', width = 22, height = 22)
# corrplot(cors, order ="AOE", type = "upper", diag = F, addCoef.col = "grey", cl.cex = 2) # 太多了看不清
corrplot::corrplot(cors, method = 'square', cl.cex = 2, tl.pos = 'lt', tl.col = 'black', tl.cex = 2.2, addCoef.col = 'grey', number.cex = 2)
dev.off()


#PCT only
pm_pct <- df1_labeled %>%
  filter(str_detect(Group, 'PCT')) %>%
  select(-(2:3)) %>%
  column_to_rownames('Filename') %>%
  t() %>%
  as.data.frame()
cors <- cor(pm_pct, use = "pairwise.complete.obs", method = 'pearson')
uptri <- cors[upper.tri(cors)]
median(uptri) # 0.9469511
min(uptri) # 0.894573
#hist(uptri)
pdf('ProteomEx_DIA_Pearson_correlation_PCT_9files_20220815.pdf', width = 22, height = 22)
# corrplot(cors, order ="AOE", type = "upper", diag = F, addCoef.col = "grey", cl.cex = 2) # 太多了看不清
corrplot::corrplot(cors, method = 'square', cl.cex = 2, tl.pos = 'lt', tl.col = 'black', tl.cex = 2.2, addCoef.col = 'grey', number.cex = 2)
dev.off()




# CV ----------------------------------------------------------------------
df_cv <- df1_labeled %>%
  group_by(Group, Amount) %>%
  summarise_at(vars(-1), function(x){sd(x, na.rm = T) / mean(x, na.rm = T)}) %>%
  ungroup() %>%
  pivot_longer(-c('Group', 'Amount'), names_to = 'prot', values_to = 'CV', values_drop_na = T)
df_cv %<>%
  add_column(Label = str_c(df_cv$Group, df_cv$Amount, sep = '-'), .before = 1)
df_cv %>% group_by(Label) %>% summarise(protnum = n()) %>%
  rio::export('ProteomEx_DIA_protein_number.xlsx')


med <- median(df_cv$CV)
p <- ggplot(df_cv, aes(x = Label, y = CV, fill = Label))+
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
ggsave('ProteomEx_DIA_CV_log2.pdf', p, width = 10, height = 5)


# raw scale
df1_labeled_raw <- df1_labeled
df1_labeled_raw[, -(1:3)] <- 2 ^ df1_labeled[, -(1:3)]

df_cv <- df1_labeled_raw %>%
  group_by(Group, Amount) %>%
  summarise_at(vars(-1), function(x){sd(x, na.rm = T) / mean(x, na.rm = T)}) %>%
  ungroup() %>%
  pivot_longer(-c('Group', 'Amount'), names_to = 'prot', values_to = 'CV', values_drop_na = T)
df_cv %<>%
  add_column(Label = str_c(df_cv$Group, df_cv$Amount, sep = '-'), .before = 1)


med <- median(df_cv$CV)
p <- ggplot(df_cv, aes(x = Label, y = CV, fill = Label))+
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
ggsave('ProteomEx_DIA_CV_raw.pdf', p, width = 10, height = 5)


