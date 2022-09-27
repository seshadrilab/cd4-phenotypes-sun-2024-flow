library(flowWorkspace)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# The purpose of this script is to compare Pneg vs. TST+ within each batch

# Create output folder
if(!dir.exists(here::here("out/QC/Batch_Effect_V2"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC/Batch_Effect_V2")))
  dir.create(here::here("out/QC/Batch_Effect_V2"), recursive = T)
}

# Load gating set
gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet"))

# Get nodes of interest (include parent nodes and markers of interest)
dput(gh_get_pop_paths(gs))
nodes <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/IL10+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD137+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD154+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CTLA4+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/OX40+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD25+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CCR7+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD39+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD73+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/IL10+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CD137+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CD154+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CTLA4+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/OX40+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/FOXP3+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CD25+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CCR7+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CD39+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CD73+")
nodes_short <- str_replace(nodes, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/", "")

# Add shortened experiment name
pData(gs)$`EXPERIMENT NAME`  <- str_replace_all(pData(gs)$`EXPERIMENT NAME`,
                                                "20220415 RSTR INS Treg B1",
                                                "B1")

pData(gs)$`EXPERIMENT NAME`  <- str_replace_all(pData(gs)$`EXPERIMENT NAME`,
                                                "20220429 RSTR INS Treg B2",
                                                "B2")

pData(gs)$`EXPERIMENT NAME` <- factor(pData(gs)$`EXPERIMENT NAME`, levels = c("B1", "B2"))

# Get counts
dmso_count <- subset(gs, Stim == "DMSO") %>%
  gs_pop_get_count_with_meta(subpopulations = nodes) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(nodes)), ~ nodes_short) 

pp1_count <- subset(gs, Stim == "PP1") %>%
  gs_pop_get_count_with_meta(subpopulations = nodes) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(nodes)), ~ nodes_short) 

tbwcl_count <- subset(gs, Stim == "TB WCL") %>%
  gs_pop_get_count_with_meta(subpopulations = nodes) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(nodes)), ~ nodes_short) 

# Plot DMSO, PP1, and TB WCL frequencies and perform Wilcoxon rank-sum test between batches
# Argument "pop" is the list of nodes of interest
# Note: p-values are unadjusted
fill_colors <- c("Pneg" = "#984EA3", "TST+" = "#4DAF4A")

plot_pop <- function(pop, counts) {     
  parent <- sub("(.*)\\/.*", "\\1", pop)
  tmp_dat <- counts %>%
    mutate(prop = !!as.name(pop) / ParentCount)
  
  ggplot(tmp_dat, aes(Status, prop)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw(base_size = 22) +
    geom_jitter(width = 0.15, height = 0, pch = 16, aes(color=!!as.symbol("Status"))) +
    labs(title = sub(".*\\/(.*)", "\\1", pop),
         subtitle = tmp_dat$Stim,
         y = sprintf("%% %s T cells", parent), 
         caption = "Treg Panel Batches") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=22),
          axis.text.y = element_text(color="black", size=15),
          axis.text.x = element_text(color="black", size=15),
          plot.title = element_text(hjust = 0.5, size=21),
          plot.subtitle = element_text(hjust = 0.5, size=13), 
          plot.caption = element_text(size=12),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(1.3, 0.2, 0, 0.2, "cm")) +
    scale_y_continuous(labels = function(x) paste0(x*100)) +
    scale_color_manual(values = fill_colors) +
    facet_wrap(~ `EXPERIMENT NAME`) +
    stat_compare_means(comparisons = list(c("Pneg", "TST+")), label = "p.format",
                       method = "wilcox.test", paired = FALSE, tip.length = 0)
}

plot_pop(nodes_short[[3]], counts = dmso_count) # The first two nodes in the list are parent nodes

for(pop in nodes_short[3:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/Batch_Effect_V2/Treg_DMSO_%s_vs_Status.png", 
                              sub("\\/", "_", pop))), width=450, height=450, units = "px")
  print(plot_pop(pop, counts = dmso_count))
  dev.off()
}

for(pop in nodes_short[3:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/Batch_Effect_V2/Treg_PP1_%s_vs_Status.png", 
                              sub("\\/", "_", pop))), width=450, height=450, units = "px")
  print(plot_pop(pop, counts = pp1_count))
  dev.off()
}

for(pop in nodes_short[3:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/Batch_Effect_V2/Treg_TBWCL_%s_vs_Status.png", 
                              sub("\\/", "_", pop))), width=450, height=450, units = "px")
  print(plot_pop(pop, counts = tbwcl_count))
  dev.off()
}
