library(here)
library(tidyverse)
library(flowWorkspace)
library(flowCore)

# The counts for the COMPASS subsets are only stored in COMPASSResult objects if they were discovered by COMPASS for that stim,
# so we have to manually add boolean gates for each subset to the GatingSets and then extract the count data later.

## Load data ##
gsPath <- here::here("out/GatingSets/RSTR_Th_GatingSet")
gs <- load_gs(gsPath)

merged_cd4_compass_data <- readRDS("processed_data/Merged_Th_CD4_COMPASS_Data.rds")
merged_cd8_compass_data <- readRDS("processed_data/Merged_Th_CD8_COMPASS_Data.rds")

## Add the CD4+ COMPASS boolean cytokine subset gates to the GatingSet ##
mapMarkers <- c("CD154", "CD137", "CTLA4", "OX40")
cd4NodeMarkerMap <- mapMarkers
# NodeMarkerMap names are gating tree paths
names(cd4NodeMarkerMap) <- paste0("CD4+", "/", c("CD154+", "CD137+", "CTLA4+", "OX40+"))

cd4_cats_mod <- as.data.frame(merged_cd4_compass_data$catsMerged) %>% 
  mutate_all(~ as.numeric(as.character(.))) %>% 
  mutate_all(~ recode(., "0" = "!", "1" = "")) %>% 
  dplyr::rename_at(vars(cd4NodeMarkerMap), ~ names(cd4NodeMarkerMap))
cd4_booleanSubsets <- cd4_cats_mod %>% 
  rowwise() %>% 
  do(booleanSubset = paste(paste0(., colnames(cd4_cats_mod)), collapse="&")) %>% 
  ungroup() %>% 
  dplyr::pull(booleanSubset) %>% 
  unlist()
names(cd4_booleanSubsets) <- paste0("CD4_", gsub("CD4\\+\\/", "", gsub("\\&", "_AND_", gsub("\\!", "NOT_", cd4_booleanSubsets))))
for(booleanSubsetName in names(cd4_booleanSubsets)) {
  # booleanSubset The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
  call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(cd4_booleanSubsets[[booleanSubsetName]])))
  g <- eval(call)
  suppressWarnings(flowWorkspace::gs_pop_add(gs, g, parent = "CD4+", name=booleanSubsetName))
}

dput(names(cd4_booleanSubsets))

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste(names(cd4_booleanSubsets), collapse="|"))))),
           parent = "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+", name = "CD4_COMPASS_Subsets")

## Add the CD8 COMPASS boolean cytokine subset gates to the GatingSet ##
mapMarkers <- c("CD154", "CD137", "CTLA4", "OX40")
cd8NodeMarkerMap <- mapMarkers
# NodeMarkerMap names are gating tree paths
names(cd8NodeMarkerMap) <- paste0("CD8+", "/", c("CD154+", "CD137+", "CTLA4+", "OX40+"))

cd8_cats_mod <- as.data.frame(merged_cd8_compass_data$catsMerged) %>%
  mutate_all(~ as.numeric(as.character(.))) %>% 
  mutate_all(~ recode(., "0" = "!", "1" = "")) %>%
  dplyr::rename_at(vars(cd8NodeMarkerMap), ~ names(cd8NodeMarkerMap))
cd8_booleanSubsets <- cd8_cats_mod %>%
  rowwise() %>%
  do(booleanSubset = paste(paste0(., colnames(cd8_cats_mod)), collapse="&")) %>%
  ungroup() %>%
  dplyr::pull(booleanSubset) %>%
  unlist()
names(cd8_booleanSubsets) <- paste0("CD8_", gsub("CD8\\+\\/", "", gsub("\\&", "_AND_", gsub("\\!", "NOT_", cd8_booleanSubsets))))
for(booleanSubsetName in names(cd8_booleanSubsets)) {
  # booleanSubset The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
  call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(cd8_booleanSubsets[[booleanSubsetName]])))
  g <- eval(call)
  suppressWarnings(flowWorkspace::gs_pop_add(gs, g, parent = "CD8+", name=booleanSubsetName))
}
dput(names(cd8_booleanSubsets))

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste(names(cd8_booleanSubsets), collapse="|"))))),
           parent = "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+", name = "CD8_COMPASS_Subsets")

## Add memory gates under CD4_COMPASS_Subsets and CD8_COMPASS_Subsets ##
cd4_ccr7_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CCR7+"
cd4_cd45ra_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD45RA+"
cd8_ccr7_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CCR7+"
cd8_cd45ra_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CD45RA+"

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = "CD4_COMPASS_Subsets", name = "Naive")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = "CD4_COMPASS_Subsets", name = "TCM")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = "CD4_COMPASS_Subsets", name = "TEMRA")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = "CD4_COMPASS_Subsets", name = "TEM")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_ccr7_path,
                                                         "&", cd8_cd45ra_path))))),
           parent = "CD8_COMPASS_Subsets", name = "Naive")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_ccr7_path,
                                                         "&!", cd8_cd45ra_path))))),
           parent = "CD8_COMPASS_Subsets", name = "TCM")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd8_ccr7_path,
                                                         "&", cd8_cd45ra_path))))),
           parent = "CD8_COMPASS_Subsets", name = "TEMRA")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd8_ccr7_path,
                                                         "&!", cd8_cd45ra_path))))),
           parent = "CD8_COMPASS_Subsets", name = "TEM")

## Add cytokine gates under CD4_COMPASS_Subsets and CD8_COMPASS_Subsets ##
cd4_ifng_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/IFNg+"
cd4_il17a_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/IL17a/IL17a+"
cd8_ifng_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/IFNg+"
cd8_il17a_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/IL17a/IL17a+"

gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd4_ifng_path),
           parent = "CD4_COMPASS_Subsets", name = "IFNg+")

gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd4_il17a_path),
           parent = "CD4_COMPASS_Subsets", name = "IL17a+")

gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd8_ifng_path),
           parent = "CD8_COMPASS_Subsets", name = "IFNg+")

gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd8_il17a_path),
           parent = "CD8_COMPASS_Subsets", name = "IL17a+")

## Add transcription factor gates under CD4_COMPASS_Subsets and CD8_COMPASS_Subsets ##
cd4_tbet_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/TBET+"
cd4_roryt_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/RORyT+"
cd8_tbet_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/TBET+"
cd8_roryt_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/RORyT+"

gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd4_tbet_path),
           parent = "CD4_COMPASS_Subsets", name = "TBET+")

gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd4_roryt_path),
           parent = "CD4_COMPASS_Subsets", name = "RORyT+")

gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd8_tbet_path),
           parent = "CD8_COMPASS_Subsets", name = "TBET+")

gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd8_roryt_path),
           parent = "CD8_COMPASS_Subsets", name = "RORyT+")

## Add chemokine gates under CD4_COMPASS_Subsets and CD8_COMPASS_Subsets ##
cd4_cxcr3_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CXCR3+"
cd4_ccr6_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CCR6+"
cd8_cxcr3_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CXCR3+"
cd8_ccr6_path <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Negative/CD8+/CCR6+"

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_cxcr3_path,
                                                         "&", cd4_ccr6_path))))),
           parent = "CD4_COMPASS_Subsets", name = "CXCR3+CCR6+")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_cxcr3_path,
                                                         "&!", cd4_ccr6_path))))),
           parent = "CD4_COMPASS_Subsets", name = "CXCR3+CCR6-")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_cxcr3_path,
                                                         "&", cd4_ccr6_path))))),
           parent = "CD4_COMPASS_Subsets", name = "CXCR3-CCR6+")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_cxcr3_path,
                                                         "&", cd8_ccr6_path))))),
           parent = "CD8_COMPASS_Subsets", name = "CXCR3+CCR6+")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_cxcr3_path,
                                                         "&!", cd8_ccr6_path))))),
           parent = "CD8_COMPASS_Subsets", name = "CXCR3+CCR6-")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd8_cxcr3_path,
                                                         "&", cd8_ccr6_path))))),
           parent = "CD8_COMPASS_Subsets", name = "CXCR3-CCR6+")

## Recompute the gating set with the new gates ##
flowWorkspace::recompute(gs)

# Plot the gating tree
png(here::here("out/Th_GatingTree_with_COMPASS_Subsets.png"), width = 7, height = 5, units = "in", res = 300)
plot(gs, bool = T, fontsize = 10)
dev.off()

# Save the gating set 
save_gs(gs, here::here("out/GatingSets/RSTR_Th_GatingSet_with_COMPASS_Subsets"))
