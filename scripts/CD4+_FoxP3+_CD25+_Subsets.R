library(tidyverse)
library(flowWorkspace)

if(!dir.exists(here::here("out/CD4+FOXP3+CD25+_Signal"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/CD4+FOXP3+CD25+_Signal")))
  dir.create(here::here("out/CD4+FOXP3+CD25+_Signal"), recursive = T)
}

# Load GatingSet
gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet/"))

# Define function
# Argument "pop" is the node of interest
plot_pop <- function(pop, freq, bg_freq, stim) {     
  parent <- sub("(.*)\\/.*", "\\1", pop)
  stim_dat <- freq %>%
    mutate(stim_prop = !!as.name(pop) / ParentCount) %>%
    select(`SAMPLE ID`, Status, stim_prop) %>%
    drop_na()
  bg_dat <- bg_freq %>%
    mutate(bg_prop = !!as.name(pop) / ParentCount) %>%
    select(`SAMPLE ID`, Status, bg_prop) %>%
    drop_na()
  bg_corr_dat <- stim_dat %>%
    left_join(bg_dat, c("SAMPLE ID", "Status")) %>%
    mutate(prop = stim_prop - bg_prop)
  p <- wilcox.test(prop ~ Status, data = bg_corr_dat, paired = FALSE)$p.value
  p.unadj.text <- sprintf("Wilcox. Rank Sum Test: p-unadj%s",
                          if_else(p < 0.001, "<0.001", paste0("=", sub("0.", ".", round(p, 3)))))
  
  ggplot(bg_corr_dat, aes(Status, prop)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw(base_size = 22) +
    geom_jitter(width = 0.15, height = 0, pch = 21, fill = "grey", alpha = 0.8) +
    labs(y = sprintf("%% %s of %s", sub(".*\\/(.*)", "\\1", pop), parent),
         caption = paste0(sprintf("%s Bg-Corr Frequencies\n", stim), p.unadj.text)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=22),
          axis.text.y = element_text(color="black", size=15),
          axis.text.x = element_text(color="black", size=15),
          plot.title = element_blank(),
          plot.caption = element_text(size=12),
          panel.grid.minor = element_blank(),
          plot.margin = margin(1.3, 0.2, 0, 0.2, "cm")) +
    scale_y_continuous(labels = function(x) paste0(x*100))
}

# Get nodes
treg_node <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+",
               "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+")
treg_node_short <- str_replace(treg_node, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/", "")

nodes <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+/IL10+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+/CD39+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+/CD73+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+/CD137+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+/CD154+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+/CCR7+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+/CTLA4+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+/OX40+")
nodes_short <- str_replace(nodes, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/CD4\\+\\/", "")

# Get Treg counts
dmso_treg_freq <- subset(gs, Stim == "DMSO") %>%
  gs_pop_get_count_with_meta(subpopulations = treg_node) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(treg_node)), ~ treg_node_short) 

pp1_treg_freq <- subset(gs, Stim == "PP1") %>%
  gs_pop_get_count_with_meta(subpopulations = treg_node) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(treg_node)), ~ treg_node_short) 

tbwcl_treg_freq <- subset(gs, Stim == "TB WCL") %>%
  gs_pop_get_count_with_meta(subpopulations = treg_node) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(treg_node)), ~ treg_node_short) 

# Get activation and functional marker counts
dmso_freq <- subset(gs, Stim == "DMSO") %>%
  gs_pop_get_count_with_meta(subpopulations = nodes) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(nodes)), ~ nodes_short) 

pp1_freq <- subset(gs, Stim == "PP1") %>%
  gs_pop_get_count_with_meta(subpopulations = nodes) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(nodes)), ~ nodes_short) 

tbwcl_freq <- subset(gs, Stim == "TB WCL") %>%
  gs_pop_get_count_with_meta(subpopulations = nodes) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(nodes)), ~ nodes_short) 

# Plot frequencies and perform Wilcoxon Rank Sum test between status groups
stim <- "PP1"
png(file=here::here(sprintf("out/CD4+FOXP3+CD25+_Signal/%s_%s_vs_Status.png", stim, 
                            sub("\\/", "_", treg_node_short[2]))), width=300, height=400, units = "px")
print(plot_pop(pop = treg_node_short[2], freq = pp1_treg_freq, bg_freq = dmso_treg_freq, stim = "PP1"))
dev.off()

stim <- "TB WCL"
png(file=here::here(sprintf("out/CD4+FOXP3+CD25+_Signal/%s_%s_vs_Status.png", stim, 
                            sub("\\/", "_", treg_node_short[2]))), width=300, height=400, units = "px")
print(plot_pop(pop = treg_node_short[2], freq = tbwcl_treg_freq, bg_freq = dmso_treg_freq, stim = "TB WCL"))
dev.off()

for(pop in nodes_short[2:length(nodes_short)]) {
  stim <- "PP1"
  png(file=here::here(sprintf("out/CD4+FOXP3+CD25+_Signal/%s_%s_vs_Status.png", stim, 
                              sub("\\/", "_", pop))), width=300, height=400, units = "px")
  print(plot_pop(pop, freq = pp1_freq, bg_freq = dmso_freq, stim = "PP1"))
  dev.off()
}

for(pop in nodes_short[2:length(nodes_short)]) {
  stim <- "TB WCL"
  png(file=here::here(sprintf("out/CD4+FOXP3+CD25+_Signal/%s_%s_vs_Status.png", stim, 
                              sub("\\/", "_", pop))), width=300, height=400, units = "px")
  print(plot_pop(pop, freq = tbwcl_freq, bg_freq = dmso_freq, stim = "TB WCL"))
  dev.off()
}
