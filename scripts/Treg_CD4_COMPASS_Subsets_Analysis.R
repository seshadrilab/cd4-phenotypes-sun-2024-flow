library(tidyverse)
library(flowWorkspace)
library(here)
library(ggplot2)
library(ggbeeswarm)
library(ggh4x)
source(here::here("scripts/Helper_Functions.R")) # for make_mag_plots() and multi_mag_plot()

# Plot background-corrected frequencies of phenotypic and functional markers among CD4 COMPASS-selected T cells

if(!dir.exists(here::here("out/Treg_CD4_COMPASS_Subsets_Analysis"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/Treg_CD4_COMPASS_Subsets_Analysis")))
  dir.create(here::here("out/Treg_CD4_COMPASS_Subsets_Analysis"), recursive = T)
}

# Load GatingSet
gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet_with_COMPASS_Subsets/"))

# Arial font setup. Downloaded afms from https://github.com/microsoft/microsoft-r-open/tree/ec3fd89e5fb5794bd8149905c134ad801bb61800
Arial <- Type1Font(family = "Arial",
                   metrics = c(here::here("data/Arial_afm/ArialMT.afm"),
                               here::here("data/Arial_afm/ArialMT-Bold.afm"),
                               here::here("data/Arial_afm/ArialMT-Italic.afm"),
                               here::here("data/Arial_afm/ArialMT-BoldItalic.afm")))
windowsFonts("Arial" = windowsFont("Arial"))
pdfFonts(Arial = Arial)

# Extract frequencies
nodes <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/FOXP3+CD25+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/FOXP3+CD25-",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/FOXP3-CD25+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/IL10+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD39+CD73+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD39+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD73+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CCR7+")
nodes_short <- str_replace(nodes, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/CD4\\+\\/CD4_COMPASS_Subsets\\/", "")

count_list <- list()

for(pop in nodes) {
  current_counts <- pData(gs) %>%
    tibble::rownames_to_column("rowname") %>%
    left_join(gs_pop_get_count_fast(gs, subpopulations = pop), by = c("rowname" = "name")) %>%
    dplyr::rename(Subpop = Count) %>%
    dplyr::select(rowname, "SAMPLE ID", "EXPERIMENT NAME", Stim, Status, Subpop, ParentCount) %>%
    mutate(Freq = (Subpop/ParentCount)*100)
  
  current_counts$Status <- factor(current_counts$Status, levels = c("Pneg", "TST+"))
  count_list[[pop]] <- current_counts
}

names(count_list) <- nodes_short

# Add column indicating "Node" in FOXP3 and CD25 count dataframes
foxp3_cd25_counts <- count_list[[nodes_short[1]]]
foxp3_cd25_counts["Node"] <- nodes_short[1]

foxp3_counts <- count_list[[nodes_short[2]]] 
foxp3_counts["Node"] <- nodes_short[2]

cd25_counts <- count_list[[nodes_short[3]]] 
cd25_counts["Node"] <- nodes_short[3]

# Combine FOXP3 and CD25 count dataframes 
treg_counts <- bind_rows(foxp3_cd25_counts, foxp3_counts) %>%
  bind_rows(cd25_counts)

nodes_compare <- c(nodes_short[1], nodes_short[2], nodes_short[3])

# Set color scheme and stims
fill_colors <- c("Pneg" = "#984EA3", "TST+" = "#4DAF4A")
stims <- c("DMSO", "PP1", "TB WCL")

# Get plots
foxp3_cd25_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[1]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[1],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 50))
                              })
names(foxp3_cd25_plots) <- stims

foxp3_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[2]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[2],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 5))
                                })
names(foxp3_plots) <- stims

cd25_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[3]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[3],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 70))
                                })
names(cd25_plots) <- stims

treg_plots <- purrr::pmap(.l = list(rep(stims, each = 2), rep(c("Pneg", "TST+"), times = 3)),
                          .f = function(n, status) {
                            multi_mag_plot(counts = treg_counts, stim_keep = n, status_keep = status, fill_colors = fill_colors, nodes_compare = nodes_compare, ylim = c(0, 45))
                            })
names(treg_plots) <- paste0(rep(stims, each = 2), "_", rep(c("Pneg", "TST+"), times = 3))

il10_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[4]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[4],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 0.3))
                                })
names(il10_plots) <- stims

cd39_cd73_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[5]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[5],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 3))
                                })
names(cd39_cd73_plots) <- stims

cd39_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[6]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[6],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 30))
                                })
names(cd39_plots) <- stims

cd73_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[7]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[7],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 35))
                                })
names(cd73_plots) <- stims

ccr7_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[8]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[8],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(70, 100))
                                })
names(ccr7_plots) <- stims

# Save plots
cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_FOXP3_CD25_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(foxp3_cd25_plots[["DMSO"]])
print(foxp3_cd25_plots[["PP1"]])
print(foxp3_cd25_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_FOXP3_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(foxp3_plots[["DMSO"]])
print(foxp3_plots[["PP1"]])
print(foxp3_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CD25_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cd25_plots[["DMSO"]])
print(cd25_plots[["PP1"]])
print(cd25_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_Treg_multi_plots.pdf"), width=6, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(treg_plots[["DMSO_Pneg"]])
print(treg_plots[["DMSO_TST+"]])
print(treg_plots[["PP1_Pneg"]])
print(treg_plots[["PP1_TST+"]])
print(treg_plots[["TB WCL_Pneg"]])
print(treg_plots[["TB WCL_TST+"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_IL10_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(il10_plots[["DMSO"]])
print(il10_plots[["PP1"]])
print(il10_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CD39_CD73_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cd39_cd73_plots[["DMSO"]])
print(cd39_cd73_plots[["PP1"]])
print(cd39_cd73_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CD39_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cd39_plots[["DMSO"]])
print(cd39_plots[["PP1"]])
print(cd39_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CD73_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cd73_plots[["DMSO"]])
print(cd73_plots[["PP1"]])
print(cd73_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CCR7_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(ccr7_plots[["DMSO"]])
print(ccr7_plots[["PP1"]])
print(ccr7_plots[["TB WCL"]])
dev.off()

# Plot and save figures
# stim <- "PP1"
# cairo_pdf(file=here::here(sprintf("out/Treg_CD4_COMPASS_Subsets_Analysis/%s_Freq_vs_Status.pdf", stim)), 
#           width=4, height=5, onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica
# for(pop in nodes_short[2:length(nodes_short)]) {
#   print(plot_pop(pop, freq = pp1_freq, bg_freq = dmso_freq, stim = stim))
#   }
# dev.off()
# 
# stim <- "TB WCL"
# cairo_pdf(file=here::here(sprintf("out/Treg_CD4_COMPASS_Subsets_Analysis/%s_Freq_vs_Status.pdf", stim)), 
#           width=4, height=5, onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica
# for(pop in nodes_short[2:length(nodes_short)]) {
#   print(plot_pop(pop, freq = tbwcl_freq, bg_freq = dmso_freq, stim = stim))
# }
# dev.off()
