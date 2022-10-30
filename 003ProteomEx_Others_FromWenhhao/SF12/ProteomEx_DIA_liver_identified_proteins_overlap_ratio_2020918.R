require(magrittr)
require(tidyverse)
library(eulerr)
library(RColorBrewer)

rm(list = ls())

# # Venn plot
# my_venn <- function(venn_list){
#   VennDiagram::venn.diagram(
#   # VennDiagram::draw.pairwise.venn(
#     x = venn_list,
#     filename = NULL,
#     resolution = 300,
#
#     compression = "lzw",
#     lwd = 1,
#     # col=c("#440154ff", '#21908dff', '#fde725ff'),
#     # fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
#     cex = 1.5,
#     fontfamily = "sans",
#     cat.cex = 1.2,
#     cat.default.pos = "outer",
#     cat.pos = c(-27, 27, 135),
#     cat.dist = c(0.055, 0.055, 0.085),
#     cat.fontfamily = "sans",
#     # cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
#     rotation = 1
#   )
# }




df <- rio::import("20220810_proteins_LR_P2P.txt") %>% column_to_rownames("prot") # log2 transformed


# remove Gel_1mm
df1 <- df %>%
  select(matches("(Gel\\d_[^1]mm)|(PCT\\d\\.\\d_\\d)")) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Filename")

# remove columns and rows with total NA
mat <- df1[, -1][, colSums(df1[, -1], na.rm = T) != 0]
df1 <- cbind(df1 %>% select(1), mat)
df1 <- df1[rowSums(df1[, -1], na.rm = T) != 0, ]

# add labels
df1_labeled <- df1 %>%
  add_column(
    Group = str_extract(df1$Filename, "(Gel)|(PCT)"),
    Amount = str_extract(df1$Filename, "(\\dmm)|(\\d\\.\\d)"),
    .before = 2
  ) %>%
  arrange(Group, Amount, Filename)

df1_labeled %<>% plyr::ddply("Amount", function(dfsub) {
  dfsub %>% add_column(
    Replicate = str_c("Replicate ", 1:nrow(dfsub)),
    .after = "Amount"
  )
})


# df_label <- df1_labeled %>% select(1:4)
# pm <- df1_labeled %>%
#   select(-(2:4)) %>%
#   column_to_rownames('Filename') %>%
#   t() %>%
#   as.data.frame()


# transform to Venn-acceptable list
df1_labeled %<>% mutate(Amount = str_c(Group, Amount, sep = ' '))
df1_venn <- df1_labeled %>% select(-Filename, -Group)
mat <- df1_venn %>% select(-(1:2))
na_position <- is.na(mat)
mat <- sapply(1:nrow(mat), function(x) {
    return(colnames(mat))
  }) %>%
  t() %>%
  as.data.frame() %>%
  setNames(colnames(mat))
mat[na_position] <- NA
df1_venn <- cbind(select(df1_venn, 1:2), mat)

venn_ls <- df1_venn %>% plyr::dlply("Amount", function(dfsub) {
  dfsub %>%
    select(-Amount) %>%
    column_to_rownames("Replicate") %>%
    t() %>%
    as.data.frame() %>%
    as.list() %>%
    lapply(function(e) e[!is.na(e)])
})

# for(i in seq_along(venn_ls)){
#   pdf(str_glue("ProteomEx_DIA_liver_identified_proteins_overlap_ratio_{names(venn_ls)[i]}_2020918.pdf"))
#   venn.plot <- my_venn(venn_ls[[i]])
#   grid::grid.draw(venn.plot)
#   graphics.off()
# }

# transform to eulerr-acceptable list
venn_overlap_ls <- lapply(venn_ls, function(e) {
  # use ratio
  universal_set_num <- e[[1]] %>%
    union(e[[2]]) %>%
    union(e[[3]]) %>%
    length()
  ratio <- c(
    "First" = e[[1]] %>% setdiff(e[[2]]) %>% setdiff(e[[3]]) %>% length(),
    "Second" = e[[2]] %>% setdiff(e[[1]]) %>% setdiff(e[[3]]) %>% length(),
    "Third" = e[[3]] %>% setdiff(e[[1]]) %>% setdiff(e[[2]]) %>% length(),
    "First&Second" = e[[1]] %>% intersect(e[[2]]) %>% setdiff(e[[3]]) %>% length(),
    "First&Third" = e[[1]] %>% intersect(e[[3]]) %>% setdiff(e[[2]]) %>% length(),
    "Second&Third" = e[[2]] %>% intersect(e[[3]]) %>% setdiff(e[[1]]) %>% length(),
    "First&Second&Third" = e[[1]] %>% intersect(e[[2]]) %>% intersect(e[[3]]) %>% length()
  ) / universal_set_num * 100
  round(ratio, 2)
})

# set color
qual_color_pals <- brewer.pal.info %>% filter(category == "qual")
all_colors <- Map(brewer.pal, qual_color_pals$maxcolors, rownames(qual_color_pals)) %>%
  unlist() %>%
  unique()
set.seed(1)
venn_colors <- sample(all_colors, 7)

# Venn Diagram
pdf("ProteomEx_DIA_liver_identified_proteins_overlap_ratio_2020919.pdf", width = 10, height = 10)
for (i in seq_along(venn_ls)) {
  p <- plot(euler(venn_overlap_ls[[i]]),
    fills = list(
      fill = c("#E0556A", "#5DF093", "#4B73E3", "#D4EB5E", "#DD75EB", "#69D1FA", "#FFFFFF"),
      edges = list(col = "white", alpha = 0),
      alpha = 0.5
    ),
    quantities = list(c(),
      col = "black",
      cex = 1.5
    ),
    labels = list(col = "black", font = 3, cex = 1.5),
    main = list(label = names(venn_overlap_ls)[i], cex = 3),
    legend = list(
      labels = names(venn_ls[[i]]),
      cex = 1.5
    )
  )
  print(p)
}
graphics.off()
