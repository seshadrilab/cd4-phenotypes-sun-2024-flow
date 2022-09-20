library(CytoML) 
library(flowCore) 
library(flowWorkspace) 
library(ggcyto)
library(here)
library(tidyverse)
library(readxl)
library(ggpubr)

## Create directories if needed ## 
if(!dir.exists(here::here("data"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC")))
  dir.create(here::here("out/QC"), recursive = T)
}

if(!dir.exists(here::here("out"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC")))
  dir.create(here::here("out/QC"), recursive = T)
}

if(!dir.exists(here::here("scripts"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC")))
  dir.create(here::here("out/QC"), recursive = T)
}

if(!dir.exists(here::here("out/QC"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC")))
  dir.create(here::here("out/QC"), recursive = T)
}

if(!dir.exists(here::here("out/QC/Counts"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC/Counts")))
  dir.create(here::here("out/QC/Counts"), recursive = T)
}

if(!dir.exists(here::here("out/QC/DMSO_Signal"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC/DMSO_Signal")))
  dir.create(here::here("out/QC/DMSO_Signal"), recursive = T)
}

if(!dir.exists(here::here("out/QC/Batch_Effect"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC/Batch_Effect")))
  dir.create(here::here("out/QC/Batch_Effect"), recursive = T)
}

if(!dir.exists(here::here("out/GatingSets"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/GatingSets")))
  dir.create(here::here("out/GatingSets"), recursive = T)
}

## Load data ##
xml_path_b1 <- here::here("data/RSTR INS Treg B1VF1.xml")
xml_path_b2 <- here::here("data/RSTR INS Treg B2VF1.xml")
fcs_subfolder_b1 <- here::here("data/20220415_RSTR_INS_Treg_FCS_B1/")
fcs_subfolder_b2 <- here::here("data/20220429_RSTR_INS_Treg_FCS_B2/")
ws_b1 <- open_flowjo_xml(xml_path_b1)
ws_b2 <- open_flowjo_xml(xml_path_b2)
metadata <- read_xlsx(here::here("data/2022 RSTR INS Metadata.xlsx"), sheet = 2, range = cell_rows(1:41))

# Drop the samples we didn't use from the metadata
metadata <- subset(metadata, is.na(metadata[,"...11"])) %>%
  select(`Sample ID`, `Status`)

## Create workspaces and prepare GatingSets ##
names(fj_ws_get_keywords(ws_b1, 117))
keywords2import <- c("EXPERIMENT NAME",
                       "$DATE",
                       "SAMPLE ID",
                       "Stim",
                       "WELL ID",
                       "PLATE NAME")

sampleGroup <- "Samples"

gs_b1 <- flowjo_to_gatingset(ws_b1,                                    
                             name=sampleGroup, 
                             keywords=keywords2import,
                             path=fcs_subfolder_b1,
                             extend_val=-10000)

gs_b2 <- flowjo_to_gatingset(ws_b2,                                    
                             name=sampleGroup, 
                             keywords=keywords2import,
                             path=fcs_subfolder_b2, 
                             extend_val=-10000)

# Check that gating trees are consistent
pop_lists_b1 <- lapply(gs_b1, gh_get_pop_paths)
unique(pop_lists_b1)

pop_lists_b2 <- lapply(gs_b2, gh_get_pop_paths)
unique(pop_lists_b2) # 3 unique gating trees!

# The second tree has an errant node at [21] and the third tree at [52]
# Which samples contain these errant nodes?
x <- as.vector(NULL)
for(i in 1:length(pop_lists_b2)){
  x[i] <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CD8+" %in% pop_lists_b2[[i]]
}
which(x == TRUE) # Entries 22-24 are true
pData(gs_b2[[c(22, 23, 24)]]) # 202060.fcs_1280948, 202062.fcs_1194574, 202064.fcs_1239513

y <- as.vector(NULL)
for(i in 1:length(pop_lists_b2)){
  y[i] <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/FSC-A, <B515-A> subset" %in% pop_lists_b2[[i]]
}
which(y == TRUE) # Entry 54 is TRUE
pData(gs_b2[[54]]) # 202130.fcs_1372548

# Remove the errant nodes from the gating set
gs_pop_remove(gs_b2[c("202060.fcs_1280948", "202062.fcs_1194574", "202064.fcs_1239513")], "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CD8+")
gs_pop_remove(gs_b2["202130.fcs_1372548"], "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/FSC-A, <B515-A> subset")

# Recheck both gating trees 
z <- c(gs_b1, gs_b2)
pop_lists <- lapply(z, gh_get_pop_paths)
unique(pop_lists)

# Remove channels from flow data that are not used by gates
gs_b1 <- gs_remove_redundant_channels(gs_b1) # drop SSC-H, G610-A, R710-A, V655-A

gs_b2 <- gs_remove_redundant_channels(gs_b2) # drop SSC-H, G610-A, R710-A, V655-A

# Add names to all channels
dput(unname(pData(parameters(gh_pop_get_data(gs_b1[[1]])))[,2]))
markernames_b1 <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CTLA4", "CD4", "OX40", "CD154", "IL10", "CD39", 
                    "FOXP3", "CD14_19", "CCR7", "CD137", "LD", "CD8a", "CD25", "CD73", "CD3")
names(markernames_b1) <- pData(parameters(gh_pop_get_data(gs_b1[[1]])))[,1]
markernames(gs_b1) <- markernames_b1
pData(parameters(gh_pop_get_data(gs_b1[[1]])))[,c(1,2)]

dput(unname(pData(parameters(gh_pop_get_data(gs_b2[[1]])))[,2]))
markernames_b2 <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CTLA4", "CD4", "OX40", "CD154", "IL10", "CD39", 
                    "FOXP3", "CD14_19", "CCR7", "CD137", "LD", "CD8a", "CD25", "CD73", "CD3")
names(markernames_b2) <- pData(parameters(gh_pop_get_data(gs_b2[[1]])))[,1]
markernames(gs_b2) <- markernames_b2
pData(parameters(gh_pop_get_data(gs_b2[[1]])))[,c(1,2)]

# Make sure nodes, pData, and markers are consistent among the two batches
setdiff(sort(gh_get_pop_paths(gs_b1)), sort(gh_get_pop_paths(gs_b2)))
all(sort(gh_get_pop_paths(gs_b1)) == sort(gh_get_pop_paths(gs_b2)))
all(markernames(gs_b1) == markernames(gs_b2))
all(colnames(pData(gs_b1)) == colnames(pData(gs_b2))) 

# Merge GatingSets from all batches
gs <- merge_list_to_gs(c(gs_b1, gs_b2))

# Add FOXP3+&CD25+ boolean gate under CD4+
cd4_foxp3_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+"
cd4_cd25_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD25+"

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_foxp3_path,
                                                         "&", cd4_cd25_path))))),
           parent = "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+", name = "FOXP3+CD25+")

# Copy gates under FOXP3+CD25+ parent gate
gates_to_copy <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/IL10+",
                   "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD39+",
                   "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD73+",
                   "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD137+",
                   "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD154+",
                   "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CCR7+",
                   "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CTLA4+",
                   "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/OX40+")

for(path in gates_to_copy) {
  gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y=path),
             parent = "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+")
}
recompute(gs)

# Plot gating tree
png(here::here("out/QC/Treg_GatingTree.png"), width = 7, height = 5, units = "in", res = 300)
plot(gs, fontsize=15, bool=T)
dev.off()

# Drop sample RS102161 (had bacterial contamination after stimulation)
gs <- subset(gs, `SAMPLE ID` != "RS102161")

# Add "Status" column indicating if a sample is Pneg or TST+
metadata <- metadata %>%
  group_by(Status) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Status, values_from = `Sample ID`) %>%
  select(-row)

pneg_samples <- na.omit(metadata$Pneg) 
pneg_samples <- paste(pneg_samples, collapse = "|")

tst_samples <- metadata$`TST+` 
tst_samples <- paste(tst_samples, collapse = "|")

pneg_index <- grepl(pneg_samples, pData(gs)$`SAMPLE ID`) 
tst_index <- grepl(tst_samples, pData(gs)$`SAMPLE ID`) 
pData(gs)$Status[pneg_index] <- "Pneg"
pData(gs)$Status[tst_index] <- "TST+"

# Save GatingSet 
save_gs(gs, here::here("out/GatingSets/RSTR_Treg_GatingSet"))

## Perform QC ##
# Load gating set if needed: 
#gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet"))

dput(gh_get_pop_paths(gs))

# Extract CD3, CD4, and CD8 counts 
cd3_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes" 
cd4_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+"
cd8_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+"
dn_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/DN"

cd3_cd4_cd8_counts <- pData(gs) %>%
  rownames_to_column(var="sample") %>%
  left_join(gs_pop_get_stats(gs, nodes = c(cd3_path, cd4_path, cd8_path, dn_path))) %>%
  pivot_wider(names_from = "pop", values_from = "count") %>%
  rename(CD3 = !!cd3_path,
         CD4 = !!cd4_path,
         CD8 = !!cd8_path,
         DN = !! dn_path)

# Plot CD3 Count
png(here::here("out/QC/Counts/CD3_Counts.png"), width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts %>% 
         mutate(Color = ifelse(cd3_path < 10000, "red", "black")),
       aes(x=`SAMPLE ID`, y=CD3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Color), position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 10000, color="red"), linetype="dashed") +
  scale_color_identity() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Plot CD4 Count
png(here::here("out/QC/Counts/CD4_Counts.png"), width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts %>% 
         mutate(Color = ifelse(cd4_path < 3000, "red", "black")),
       aes(x=`SAMPLE ID`, y=CD4)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Color), position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 3000, color="red"), linetype="dashed") +
  scale_color_identity() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Plot CD8 Count
png(here::here("out/QC/Counts/CD8_Counts.png"), width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts %>% 
         mutate(Color = ifelse(cd8_path < 3000, "red", "black")),
       aes(x=`SAMPLE ID`, y=CD8)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Color), position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 3000, color="red"), linetype="dashed") +
  scale_color_identity() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Plot DN Count
png(here::here("out/QC/Counts/DN_Counts.png"), width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts %>% 
         mutate(Color = ifelse(dn_path < 3000, "red", "black")),
       aes(x=`SAMPLE ID`, y=DN)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Color), position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 3000, color="red"), linetype="dashed") +
  scale_color_identity() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Investigate the low CD4 and CD8 count wells
low_count <- cd3_cd4_cd8_counts %>%
  dplyr::filter(CD4 < 3000 | CD8 < 3000 | DN < 3000) %>%
  select("SAMPLE ID", "Stim", "CD4", "CD8", "DN") %>%
  arrange("SAMPLE ID")

low_count # None for CD4 and CD8, but 7 samples have low DN count. Nothing to drop for COMPASS

## Plot DMSO signal stratified by status ##
# Load gating set if needed: 
# gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet"))

# Get nodes of interest (include parent nodes and markers of interest)
dput(gh_get_pop_paths(gs))
nodes <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/DN",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/IL10+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD137+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD154+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CTLA4+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/OX40+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/IL10+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CD137+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CD154+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/CTLA4+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD8+/OX40+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/DN/IL10+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/DN/CD137+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/DN/CD154+", 
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/DN/CTLA4+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/DN/OX40+")
nodes_short <- str_replace(nodes, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/", "")

# Get DMSO counts
dmso_count <- subset(gs, Stim == "DMSO") %>%
  gs_pop_get_count_with_meta(subpopulations = nodes) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(nodes)), ~ nodes_short) 

# Plot DMSO frequencies and perform Wilcoxon rank-sum test between status groups
# Argument "pop" is the list of nodes of interest
plot_pop <- function(pop) {     
  parent <- sub("(.*)\\/.*", "\\1", pop)
  tmp_dat <- dmso_count %>%
    mutate(prop = !!as.name(pop) / ParentCount)
  wilcox_p <- wilcox.test(prop ~ Status, data = tmp_dat, paired = FALSE)$p.value
  p.unadj.text <- sprintf("Wilcoxon Rank-Sum Test: p-unadj%s",
                          if_else(wilcox_p < 0.001, "<0.001", paste0("=", sub("0.", ".", round(wilcox_p, 3)))))
  
  ggplot(tmp_dat, aes(Status, prop)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw(base_size = 22) +
    geom_jitter(width = 0.15, height = 0, pch = 21, fill = "grey", alpha = 0.8) +
    labs(y = sprintf("%% %s of %s", sub(".*\\/(.*)", "\\1", pop), parent),
         caption = paste0("DMSO frequencies\n", p.unadj.text)) +
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

plot_pop(nodes_short[[4]]) # The first three nodes in the list are parent nodes

for(pop in nodes_short[4:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/DMSO_Signal/Treg_%s_vs_Status.png", 
                              sub("\\/", "_", pop))), width=300, height=265, units = "px")
  print(plot_pop(pop))
  dev.off()
}

## Check for batch effect ##
# Load gating set if needed: 
# gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet"))

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
  # wilcox_p <- wilcox.test(prop ~ `EXPERIMENT NAME`, data = tmp_dat, paired = FALSE)$p.value
  # p.unadj.text <- sprintf("Wilcoxon Rank-Sum Test: p-unadj%s",
  #                         if_else(wilcox_p < 0.001, "<0.001", paste0("=", sub("0.", ".", round(wilcox_p, 3)))))
  
  ggplot(tmp_dat, aes(`EXPERIMENT NAME`, prop)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw(base_size = 22) +
    geom_jitter(width = 0.15, height = 0, pch = 16, aes(color=!!as.symbol("Status"))) +
    labs(title = tmp_dat$Stim,
         y = sprintf("%% %s of %s", sub(".*\\/(.*)", "\\1", pop), parent), 
         caption = "Treg Panel Batches") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=22),
          axis.text.y = element_text(color="black", size=15),
          axis.text.x = element_text(color="black", size=15),
          plot.title = element_text(size=17, hjust=0.5),
          plot.caption = element_text(size=12),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(1.3, 0.2, 0, 0.2, "cm")) +
    scale_y_continuous(labels = function(x) paste0(x*100)) +
    scale_color_manual(values = fill_colors) +
    facet_wrap(~ Status) +
    stat_compare_means(comparisons = list(c("B1", "B2")), label = "p.format",
                       method = "wilcox.test", paired = FALSE, tip.length = 0)
}

plot_pop(nodes_short[[3]], counts = dmso_count) # The first two nodes in the list are parent nodes

for(pop in nodes_short[3:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/Batch_Effect/Treg_DMSO_%s_vs_Batch.png", 
                              sub("\\/", "_", pop))), width=450, height=450, units = "px")
  print(plot_pop(pop, counts = dmso_count))
  dev.off()
}

for(pop in nodes_short[3:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/Batch_Effect/Treg_PP1_%s_vs_Batch.png", 
                              sub("\\/", "_", pop))), width=450, height=450, units = "px")
  print(plot_pop(pop, counts = pp1_count))
  dev.off()
}

for(pop in nodes_short[3:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/Batch_Effect/Treg_TBWCL_%s_vs_Batch.png", 
                              sub("\\/", "_", pop))), width=450, height=450, units = "px")
  print(plot_pop(pop, counts = tbwcl_count))
  dev.off()
}
