library(tidyverse)
library(magrittr)
require(ggsignif)
require(ggpubr)

rm(list = ls())

# c('#064F89', '#E6B429', '#8FC31F', '#DD704A')
my_colors <- c(`In-solution` = '#064F89', 
               PCT = '#E6B429',
               `proExM-MS` = '#8FC31F', 
               ProteomEx = '#DD704A', `ProteomEx-2` = '#DD704A'
               )


# barplot function
my_barplot <- function(df_rlt){
  df_bar <- df_rlt %>%
    group_by(Group) %>%
    summarise(
      average = mean(Ratio, na.rm = T),
      sd = sd(Ratio, na.rm = T)
    )
  
  this_colors <- my_colors[sapply(names(my_colors), function(e) e %in% df_bar$Group)]
  # barplot
  p1 <- ggplot(data = df_bar) +
    geom_bar(aes(x = Group, color = Group, weight = average), fill = 'white') +
    geom_errorbar(aes(x = Group, y = average, ymin = average - sd, ymax = average + sd, color = Group), width = 0.5)+
    # geom_point(aes(x = Group, y = Ratio), data = df_rlt, shape = "circle", size = 3) +
    # geom_jitter(aes(x = Group, y = Ratio), data = df_rlt)+
    # scale_fill_brewer(palette = "Set2") +
    # scale_color_brewer(palette = "Set2") +
    scale_fill_manual(values = this_colors) +
    scale_color_manual(values = this_colors) +
    labs(x = "",
         y = "Anchor modified peptide ratio (%)") +
    # theme_minimal()+
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.line = element_line(colour = "black"),
      panel.background = element_blank(),
    )
  
  if('ProteomEx-2' %in% df_bar$Group){
    p1 <- p1 + geom_signif(aes(x = Group, y = Ratio),
                           data = df_rlt,
                           comparisons = list(
                             c('In-solution', 'ProteomEx'),
                             c('PCT', 'ProteomEx'),
                             c('proExM-MS', 'ProteomEx'),
                             c('In-solution', 'ProteomEx-2'),
                             c('PCT', 'ProteomEx-2'),
                             c('proExM-MS', 'ProteomEx-2'),
                             c('ProteomEx', 'ProteomEx-2')
                           ),
                           step_increase = 0.1,
                           #map_signif_level = T,
                           # test = function(x1, x2) t.test(x1, x2, alternative = 'less')) # 单尾
                           test = t.test) # 双尾,
  } else {
    p1 <- p1 + geom_signif(aes(x = Group, y = Ratio),
                           data = df_rlt,
                           comparisons = list(
                             c('In-solution', 'ProteomEx'),
                             c('PCT', 'ProteomEx'),
                             c('proExM-MS', 'ProteomEx')
                           ),
                           step_increase = 0.1,
                           #map_signif_level = T,
                           # test = function(x1, x2) t.test(x1, x2, alternative = 'less')) # 单尾
                           test = t.test) # 双尾
  }
  
  return(p1)
}



# main --------------------------------------------------------------------

quantbl_modrat <- function(df, pattern){
  ##-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # certain modification ratio of peptide sequences in a quant_csv table
  # df: table read from FragPipe IonQuant `_quant.csv`
  # pattern: modification such like [54.0474]; 
  #          should be in regular expression format
  ##-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # statistic
  x1 <- df %>% nrow
  x2 <- df %>% filter(str_detect(`Assigned Modifications`, pattern)) %>% nrow
  ratio <- x2 / x1 * 100 # %

  return(ratio)
}

# read tables
tsvs <- list.files('//172.16.13.136/ProteomEx/ProteomEx_LFQ_20220606/para_search/ProteomEx_3variable_LFQ/',
                          pattern = '^peptide\\.tsv$',
                          full.names = T,
                          recursive = T)

# tsvs <- list.files('//172.16.13.136/ProteomEx/ProteomEx_LFQ_20220606/para_search/AC_LFQ',
#                       pattern = 'peptide\\.tsv$',
#                       full.names = T,
#                       recursive = T) # AC mod
# 
# 
# # get modification ratio
promex_ratio <- sapply(tsvs, function(tsv){
  df <- read.delim(tsv, check.names = F)
  ratio <- quantbl_modrat(df, '\\(54.0474\\)')
  return(ratio)
})
names(promex_ratio) <- str_split(tsvs, '/') %>% lapply(tail, 2) %>% sapply(head, 1) %>% str_replace('_T$', '')

ac_ratio <- sapply(tsvs, function(tsv){
  df <- read.delim(tsv, check.names = F)
  ratio <- quantbl_modrat(df, '\\(114.1656\\)|\\(168.2130\\)')
  return(ratio)
})
ac_ratio1 <- sapply(tsvs, function(tsv){
  df <- read.delim(tsv, check.names = F)
  ratio <- quantbl_modrat(df, '\\(114.1656\\)')
  return(ratio)
})
ac_ratio2 <- sapply(tsvs, function(tsv){
  df <- read.delim(tsv, check.names = F)
  ratio <- quantbl_modrat(df, '\\(168.2130\\)')
  return(ratio)
})
names(ac_ratio) <- names(ac_ratio1) <- names(ac_ratio2) <- str_split(tsvs, '/') %>% lapply(tail, 2) %>% sapply(head, 1) %>% str_replace('_AC$', '')



df_rlt <- data.frame(
  File = names(promex_ratio),
  `Ratio X` = promex_ratio,
  `Ratio Y` = ac_ratio,
  `Ratio Y1` = ac_ratio1,
  `Ratio Y2` = ac_ratio2,
  check.names = F
)
df_rlt %<>% add_column(Group = str_extract(df_rlt$File, '[abcgS]\\d') %>%
                   str_to_upper() %>%
                   str_replace_all('_', '') %>%
                   str_replace_all('\\d', ''), #%>%
                   # str_replace('^S$', 'A'), # S should be A
                 .after = 'File')

# str(df_rlt)

df_rlt$Group %<>% sapply(function(e){
  switch(e, 
         A = 'ProteomEx',
         S = 'ProteomEx-2',
         B = 'PCT',
         C = 'In-solution',
         G = 'proExM-MS')
}) %>% unname()

df_rlt$Group %<>% factor(levels = c('In-solution', 'PCT', 'proExM-MS', 'ProteomEx', 'ProteomEx-2'),
                         ordered = T)

df_rlt %<>% arrange(Group) %>%
  pivot_longer(cols = -(1:2), names_to = 'Modification', values_to = 'Ratio') %>%
  mutate(Modification = str_replace(Modification, 'Ratio ', ''))
df_rlt$Modification %<>% sapply(function(e){
  switch(e, 
         X = '54.0474',
         Y = '114.1656 + 168.2130',
         Y1 = '114.1656',
         Y2 = '168.2130')
}) %>% unname()
df_rlt$Modification %<>% factor(levels = c('54.0474', '114.1656', '168.2130', '114.1656 + 168.2130'),
                         ordered = T)

# p1 <- my_barplot(df_rlt)
df_bar <- df_rlt %>%
  group_by(Group, Modification) %>%
  summarise_if(is.numeric, list(average = mean, sd = sd)) %>%
  ungroup() %>%
  left_join(df_rlt, ., by = c('Group', 'Modification')) %>%
  arrange(Modification, Group)

this_colors <- my_colors[sapply(names(my_colors), function(e) e %in% df_bar$Group)]

# change df_bar for plotting
gmpair_ls <- df_bar %>% distinct(Modification, Group) %>%
  t() %>%
  as.data.frame() %>%
  as.list()
for(i in seq_along(gmpair_ls)){
  gmpair <- gmpair_ls[[i]]
  matched_pos <- which(df_bar$Group == gmpair[1] & df_bar$Modification == gmpair[2])
  for(j in matched_pos[-1]){
    df_bar[j, 'average'] <- NaN
    df_bar[j, 'sd'] <- NaN
  }
}


# barplot
p <- ggplot(data = df_bar) +
  facet_grid(. ~Modification) +
  geom_bar(aes(x = Group, color = Group, weight = average),
           data = df_bar %>% select(-File, Ratio) %>% distinct(),
           position = 'dodge', fill = 'white') +
  geom_text(
    aes(label = round(average, 2), x = Group, y = max * 1.02),
    data = df_bar %>% group_by(Modification, Group) %>% summarise(max = max(Ratio), .groups = 'drop') %>% left_join(df_bar) %>% drop_na(),
    position = position_dodge(0.9), size = 4,
    vjust = 0
  )+
  geom_errorbar(aes(x = Group, y = average,
                    ymin = average - sd, ymax = average + sd, color = Group),
                data = distinct(df_bar, Group, Modification, average, sd),
                width = 0.5, position = 'dodge')+
  # geom_point(aes(x = Group, y = Ratio), shape = "circle", size = 3) +
  geom_jitter(aes(x = Group, y = Ratio, color = Group), size = 2)+
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = this_colors) +
  scale_color_manual(values = this_colors) +
  labs(x = "",
       y = "Anchor modified peptide ratio (%)") +
  # theme_minimal()+
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
  )



p1 <- p + 
  geom_signif(aes(x = Group, y = Ratio),
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
              test = t.test) # 双尾

ggsave('ProteomEx_LFQ_anchor_modification_ratio_SF6A_20221028.pdf', p1, width = 16, height = 9)

df_bar %>%
  arrange(Modification, Group) %>%
  rio::export('ProteomEx_LFQ_anchor_modification_ratio_SF6A_20221028.xlsx')


p2 <- p + 
  geom_signif(aes(x = Group, y = Ratio),
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

ggsave('ProteomEx_LFQ_anchor_modification_ratio_SF6A_signifmark_20221028.pdf', p2, width = 16, height = 9)
