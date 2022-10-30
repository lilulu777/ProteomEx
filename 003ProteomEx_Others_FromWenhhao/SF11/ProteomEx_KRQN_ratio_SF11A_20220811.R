# use variable modifications of both AC and ProteomEx
library(tidyverse)
library(magrittr)
# require(ggsignif)

rm(list = ls())

my_colors <- c(`In-solution` = '#064F89', 
               PCT = '#E6B429',
               `proExM-MS` = '#8FC31F', 
               ProteomEx = '#DD704A',
               `ProteomEx-2` = '#DD704A')

group_rename <- c('In-solution', 'PCT', 'proExM-MS', 'ProteomEx')
names(group_rename) <- c('c', 'b', 'g', 'a')


files <- list.files('//172.16.13.136/ProteomEx/ProteomEx_LFQ_20220606/para_search/ProteomEx_3variable_LFQ/', pattern = '^peptide\\.tsv$', full.names = T, recursive = T)

rlt_ls <- lapply(files, function(f){
  df <- read.delim(f, check.names = F)
  data.frame(
    Ratio = sapply(c('K', 'R', 'N', 'Q'), function(pattern){
      sum(str_detect(df$Peptide, pattern))
    }) / nrow(df)
  )
})

df_rlt <- Reduce(cbind, rlt_ls)
colnames(df_rlt) <- str_split(files, '/') %>% lapply(tail, 2) %>% sapply(head, 1) %>% str_replace('_T$', '')

df_rlt %<>%
  `*`(100) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('File') %>%
  pivot_longer(-File, names_to = 'Amino acid', values_to = 'Peptide ratio (%)') %>%
  mutate(Group = str_extract(File, '[abcSg]\\d') %>% str_replace_all('_', '') %>% str_replace_all('\\d', ''))


df_rlt$Group %<>% sapply(function(e){
  switch(e, 
         a = 'ProteomEx',
         S = 'ProteomEx-2',
         b = 'PCT',
         c = 'In-solution',
         g = 'proExM-MS')}) %>% unname()

df_rlt$Group %<>% factor(levels = c('In-solution', 'PCT', 'proExM-MS', 'ProteomEx', 'ProteomEx-2'),
                         ordered = T)
df_rlt %<>% arrange(Group)
df_bar <- df_rlt %>%
  group_by(Group, `Amino acid`) %>%
  summarise(
    average = mean(`Peptide ratio (%)`, na.rm = T),
    sd = sd(`Peptide ratio (%)`, na.rm = T),
    .groups = 'drop'
  ) %>%
  left_join(df_rlt, ., by = c('Group', 'Amino acid')) %>%
  arrange(`Amino acid`, Group)


this_colors <- my_colors[sapply(names(my_colors), function(e) e %in% df_bar$Group)]

gmpair_ls <- df_bar %>% distinct(`Amino acid`, Group) %>%
  t() %>%
  as.data.frame() %>%
  as.list()
for(i in seq_along(gmpair_ls)){
  gmpair <- gmpair_ls[[i]]
  matched_pos <- which(df_bar$Group == gmpair[1] & df_bar$`Amino acid` == gmpair[2])
  for(j in matched_pos[-1]){
    df_bar[j, 'average'] <- NaN
    df_bar[j, 'sd'] <- NaN
  }
}


p <- ggplot(data = df_bar) +
  facet_grid(. ~`Amino acid`) +
  geom_bar(aes(x = Group, color = Group, weight = average), fill = 'white', position = 'dodge') +
  geom_text(
    aes(label = round(average, 2), x = Group, y = max * 1.02),
    data = df_bar %>% group_by(`Amino acid`, Group) %>% summarise(max = max(`Peptide ratio (%)`), .groups = 'drop') %>% left_join(df_bar) %>% drop_na(),
    position = position_dodge(0.9), size = 2,
    vjust = 0
  )+
  geom_errorbar(aes(x = Group, y = average, ymin = average - sd, ymax = average + sd, color = Group), position = position_dodge(0.9), width = 0.5)+
  geom_jitter(aes(x = Group, y = `Peptide ratio (%)`, color = Group), size = 2)+
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = this_colors) +
  scale_color_manual(values = this_colors) +
  labs(x = "Amino acid", y = "Peptide with certain amino acid ratio (%)") +
  # theme_minimal()+
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
  )


p1 <- p + 
  geom_signif(aes(x = Group, y = `Peptide ratio (%)`),
              data = df_rlt,
              comparisons = combn(unique(df_bar$Group), 2, simplify = F),
              # comparisons = list(
              #   c('In-solution', 'ProteomEx'),
              #   c('PCT', 'ProteomEx'),
              #   c('proExM-MS', 'ProteomEx')
              # ),
              step_increase = 0.1,
              #map_signif_level = T,
              map_signif_level = function(p) sprintf("%.3g", p),
              # test = function(x1, x2) t.test(x1, x2, alternative = 'less')) # 单尾
              test = t.test)
ggsave('ProteomEx_KRQN_ratio_SF11A_20221028.pdf', p1, width = 8, height = 4)
ggsave('ProteomEx_KRQN_ratio_SF11A_20221028_view.pdf', p1, width = 16, height = 16)

df_bar %>%
  arrange(`Peptide ratio (%)`, Group) %>%
  rio::export('ProteomEx_KRQN_ratio_SF11A_20221028.xlsx')



p2 <- p + 
  geom_signif(aes(x = Group, y = `Peptide ratio (%)`),
              data = df_rlt,
              comparisons = combn(unique(df_bar$Group), 2, simplify = F),
              # comparisons = list(
              #   c('In-solution', 'ProteomEx'),
              #   c('PCT', 'ProteomEx'),
              #   c('proExM-MS', 'ProteomEx')
              # ),
              step_increase = 0.1,
              map_signif_level = function(p) {
                if(p < 0.001) '***'
                else if(p < 0.01) '**'
                else if (p < 0.05) '*'
                else ''
              },
              # test = function(x1, x2) t.test(x1, x2, alternative = 'less')) # 单尾
              test = t.test) # 双尾
ggsave('ProteomEx_KRQN_ratio_SF11A_signifmark_20221028_view.pdf', p2, width = 16, height = 16)
ggsave('ProteomEx_KRQN_ratio_SF11A_signifmark_20221028.pdf', p2, width = 8, height = 4)




# only a/b/c/g
df_rlt %<>% filter(!(Group %in% c('ProteomEx-2')))
df_bar <- df_rlt %>%
  group_by(Group, `Amino acid`) %>%
  summarise(
    average = mean(`Peptide ratio (%)`, na.rm = T),
    sd = sd(`Peptide ratio (%)`, na.rm = T),
    .groups = 'drop'
  ) %>%
  left_join(df_rlt, ., by = c('Group', 'Amino acid')) %>%
  arrange(`Amino acid`, Group)


this_colors <- my_colors[sapply(names(my_colors), function(e) e %in% df_bar$Group)]

gmpair_ls <- df_bar %>% distinct(`Amino acid`, Group) %>%
  t() %>%
  as.data.frame() %>%
  as.list()
for(i in seq_along(gmpair_ls)){
  gmpair <- gmpair_ls[[i]]
  matched_pos <- which(df_bar$Group == gmpair[2] & df_bar$`Amino acid` == gmpair[1])
  for(j in matched_pos[-1]){
    df_bar[j, 'average'] <- NaN
    df_bar[j, 'sd'] <- NaN
  }
}


p <- ggplot(data = df_bar) +
  facet_grid(. ~`Amino acid`) +
  geom_bar(aes(x = Group, color = Group, weight = average), fill = 'white', position = 'dodge') +
  geom_text(
    aes(label = round(average, 2), x = Group, y = max * 1.02),
    data = df_bar %>% group_by(`Amino acid`, Group) %>% summarise(max = max(`Peptide ratio (%)`), .groups = 'drop') %>% left_join(df_bar) %>% drop_na(),
    position = position_dodge(0.9), size = 2,
    vjust = 0
  )+
  geom_errorbar(aes(x = Group, y = average, ymin = average - sd, ymax = average + sd, color = Group), position = position_dodge(0.9), width = 0.5)+
  geom_jitter(aes(x = Group, y = `Peptide ratio (%)`, color = Group), size = 2)+
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = this_colors) +
  scale_color_manual(values = this_colors) +
  labs(x = "Amino acid", y = "Peptide with certain amino acid ratio (%)") +
  # theme_minimal()+
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
  )


p1 <- p + 
  geom_signif(aes(x = Group, y = `Peptide ratio (%)`),
              data = df_rlt,
              comparisons = combn(unique(df_bar$Group), 2, simplify = F),
              # comparisons = list(
              #   c('In-solution', 'ProteomEx'),
              #   c('PCT', 'ProteomEx'),
              #   c('proExM-MS', 'ProteomEx')
              # ),
              step_increase = 0.1,
              #map_signif_level = T,
              map_signif_level = function(p) sprintf("%.3g", p),
              # test = function(x1, x2) t.test(x1, x2, alternative = 'less')) # 单尾
              test = t.test)
ggsave('ProteomEx_KRQN_ratio_SF11A_abcgOnly_20221028.pdf', p1, width = 8, height = 4)
ggsave('ProteomEx_KRQN_ratio_SF11A_abcgOnly_20221028_view.pdf', p1, width = 16, height = 16)

df_bar %>%
  arrange(`Peptide ratio (%)`, Group) %>%
  rio::export('ProteomEx_KRQN_ratio_SF11A_abcgOnly_20221028.xlsx')



p2 <- p + 
  geom_signif(aes(x = Group, y = `Peptide ratio (%)`),
              data = df_rlt,
              comparisons = combn(unique(df_bar$Group), 2, simplify = F),
              # comparisons = list(
              #   c('In-solution', 'ProteomEx'),
              #   c('PCT', 'ProteomEx'),
              #   c('proExM-MS', 'ProteomEx')
              # ),
              step_increase = 0.1,
              map_signif_level = function(p) {
                if(p < 0.001) '***'
                else if(p < 0.01) '**'
                else if (p < 0.05) '*'
                else ''
              },
              # test = function(x1, x2) t.test(x1, x2, alternative = 'less')) # 单尾
              test = t.test) # 双尾
ggsave('ProteomEx_KRQN_ratio_SF11A_signifmark_abcgOnly_20221028_view.pdf', p2, width = 16, height = 16)
ggsave('ProteomEx_KRQN_ratio_SF11A_signifmark_abcgOnly_20221028.pdf', p2, width = 8, height = 4)



df_bar %<>% arrange(`Amino acid`, Group)

# t-test
an_vec <- unique(df_rlt$`Amino acid`)
pvalue_ls <- list()
for(i in seq_along(an_vec)){
  an <- an_vec[i]
  pair_ls <- combn(unique(df_rlt$Group), 2, simplify = F)
  rlt <- sapply(pair_ls, function(pair){
    dfsub <- df_rlt %>% filter(`Amino acid` == an, Group %in% pair)
    x1 <- dfsub %>% filter(Group == pair[1]) %>% pull(`Peptide ratio (%)`)
    x2 <- dfsub %>% filter(Group == pair[2]) %>% pull(`Peptide ratio (%)`)
    t.test(x1, x2, var.equal = F)$p.value
  })
  names(rlt) <- sapply(pair_ls, function(e) str_c(e, collapse = ' versus '))
  pvalue_ls[[i]] <- rlt
  
}
df_p <- Reduce(cbind, pvalue_ls)
colnames(df_p) <- an_vec
df_p2 <- df_p %>% as.data.frame() %>% rownames_to_column('T-test')
df_p1 <- df_p2 %>% pivot_longer(-`T-test`, names_to = 'Amino acid', values_to = 'p-value')
rio::export(list(df_p1, df_p2), 'ProteomEx_peptide_ratio_with_certain_amino_acid_20220811.xlsx')


list(stat = df_bar,
     p_value1 = df_p1,
     p_value2 = df_p2) %>% rio::export('ProteomEx_KRQN_ratio_SF11_v2_20220811.xlsx')
