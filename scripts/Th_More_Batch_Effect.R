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
gs <- load_gs(here::here("out/GatingSets/RSTR_Th_GatingSet"))

# Add IL17a+ gate directly under CD4+ and CD8+
cd4_il17a_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/IL17a/IL17a+"
cd8_il17a_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/IL17a/IL17a+"

gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd4_il17a_path),
           parent = "CD4+", name = "IL17a+")
gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd8_il17a_path),
           parent = "CD8+", name = "IL17a+")

recompute(gs)

# Get nodes of interest (include parent nodes and markers of interest)
dput(gh_get_pop_paths(gs))
nodes <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD45RA+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CCR7+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CXCR3+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CCR6+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/TBET+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/RORyT+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/IFNg+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/IL17a+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD137+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD154+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CTLA4+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/OX40+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CD45RA+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CCR7+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CXCR3+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CCR6+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/TBET+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/RORyT+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/IFNg+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/IL17a+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CD137+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CD154+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CTLA4+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/OX40+")
nodes_short <- str_replace(nodes, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/CD4 Positive\\/", "")
nodes_short <- str_replace(nodes_short, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/CD4 Negative\\/", "")

# Add shortened experiment name
pData(gs)$`EXPERIMENT NAME`  <- str_replace_all(pData(gs)$`EXPERIMENT NAME`,
                                                "20220415 RSTR INS Th B1",
                                                "B1")

pData(gs)$`EXPERIMENT NAME`  <- str_replace_all(pData(gs)$`EXPERIMENT NAME`,
                                                "20220429 RSTR INS Th B2",
                                                "B2")

pData(gs)$Batch <- factor(pData(gs)$`EXPERIMENT NAME`, levels = c("B1", "B2"))

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
         caption = "Th Panel Batches") +
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
    facet_wrap(~ Batch) +
    stat_compare_means(comparisons = list(c("Pneg", "TST+")), label = "p.format",
                       method = "wilcox.test", paired = FALSE, tip.length = 0)
}

plot_pop(nodes_short[[3]], counts = dmso_count) # The first two nodes in the list are parent nodes

for(pop in nodes_short[3:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/Batch_Effect_V2/Th_DMSO_%s_vs_Status.png", 
                              sub("\\/", "_", pop))), width=450, height=450, units = "px")
  print(plot_pop(pop, counts = dmso_count))
  dev.off()
}

for(pop in nodes_short[3:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/Batch_Effect_V2/Th_PP1_%s_vs_Status.png", 
                              sub("\\/", "_", pop))), width=450, height=450, units = "px")
  print(plot_pop(pop, counts = pp1_count))
  dev.off()
}

for(pop in nodes_short[3:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/Batch_Effect_V2/Th_TBWCL_%s_vs_Status.png", 
                              sub("\\/", "_", pop))), width=450, height=450, units = "px")
  print(plot_pop(pop, counts = tbwcl_count))
  dev.off()
}
