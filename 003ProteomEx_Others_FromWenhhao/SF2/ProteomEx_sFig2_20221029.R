pacman::p_unload(pacman::p_loaded(), character.only = T)
library(tidyverse)
library(magrittr)
require(ggsignif)

rm(list = ls())


df <- rio::import('ProteomEx_sFig2_input_20221029.xlsx')
df %<>% mutate(Group = str_c(`homogenization buffer`, hydrogel, `chemical anchors`, sep = '+')) %>%
  select(-hydrogel, -`chemical anchors`, -`homogenization buffer`)

df_rlt_reg <- df_rlt <- df %>%
  group_by(Group) %>%
  summarise(average_pep = mean(peptides),
            sd_pep = sd(peptides),
            average_prot = mean(proteins),
            sd_prot = sd(proteins)) %>%
  left_join(df, ., by = c('Group'))
df_rlt$Group %<>% factor(levels = c('SDS+PAE+NSA','TFE+PAE+NSA','SDS+TPT+AGE','TFE+TPT+AGE','SDS+PAE+NAS','TFE+PAE+NAS'), ordered = T)


# change data for plotting
groups <- levels(df_rlt$Group)
for(i in seq_along(groups)){
  grp <- groups[i]
  matched_pos <- which(df_rlt$Group == grp)
  for(j in matched_pos[-1]){
    df_rlt[j, 'average_pep'] <- NaN
    df_rlt[j, 'sd_pep'] <- NaN
    df_rlt[j, 'average_prot'] <- NaN
    df_rlt[j, 'sd_prot'] <- NaN
  }
}

this_colors <- c('#064F89', '#064F8955', '#E6B429', '#E6B42955', '#8FC31F', '#8FC31F55')


# 1.peptide ---------------------------------------------------------------
p1 <- ggplot(data = df_rlt) +
  geom_bar(aes(x = Group, color = Group, weight = average_pep), fill = 'white', position = 'dodge') +
  geom_text(
    aes(label = round(average_pep, 0), x = Group, y = max * 1.02),
    data = df_rlt %>% group_by(Group) %>% summarise(max = max(peptides)) %>% left_join(df_rlt, by = 'Group') %>% drop_na(),
    position = position_dodge(0.9), size = 4,
    vjust = 0
  )+
  geom_errorbar(aes(x = Group, y = average_pep, ymin = average_pep - sd_pep, ymax = average_pep + sd_pep, color = Group), position = position_dodge(0.9), width = 0.5)+
  geom_jitter(aes(x = Group, y = peptides, color = Group), size = 2)+
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = this_colors) +
  scale_color_manual(values = this_colors) +
  labs(x = "", y = "Peptides") +
  # theme_minimal()+
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
  )
ggsave('ProteomEx_SF2_optimization_of_tissue_expansion_protocol_peptide.pdf', p1, width = 8, height = 4)



# 2.protein ---------------------------------------------------------------
p2 <- ggplot(data = df_rlt) +
  geom_bar(aes(x = Group, color = Group, weight = average_prot), fill = 'white', position = 'dodge') +
  geom_text(
    aes(label = round(average_prot, 0), x = Group, y = max * 1.02),
    data = df_rlt %>% group_by(Group) %>% summarise(max = max(proteins)) %>% left_join(df_rlt, by = 'Group') %>% drop_na(),
    position = position_dodge(0.9), size = 4,
    vjust = 0
  )+
  geom_errorbar(aes(x = Group, y = average_prot, ymin = average_prot - sd_prot, ymax = average_prot + sd_prot, color = Group), position = position_dodge(0.9), width = 0.5)+
  geom_jitter(aes(x = Group, y = proteins, color = Group), size = 2)+
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = this_colors) +
  scale_color_manual(values = this_colors) +
  labs(x = "", y = "Proteins") +
  # theme_minimal()+
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
  )
ggsave('ProteomEx_SF2_optimization_of_tissue_expansion_protocol_protein.pdf', p2, width = 8, height = 4)

