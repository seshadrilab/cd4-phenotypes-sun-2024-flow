library(tidyverse)
library(flowWorkspace)
library(here)
library(ggplot2)
library(ggbeeswarm)
library(ggh4x)
source(here::here("scripts/Helper_Functions.R")) # for make_mag_plots() and multi_mag_plot()

# Plot background-corrected frequencies of phenotypic and functional markers among CD4 COMPASS-selected T cells

if(!dir.exists(here::here("out/Th_CD4_COMPASS_Subsets_Analysis"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/Th_CD4_COMPASS_Subsets_Analysis")))
  dir.create(here::here("out/Th_CD4_COMPASS_Subsets_Analysis"), recursive = T)
}

# Load GatingSet
gs <- load_gs(here::here("out/GatingSets/RSTR_Th_GatingSet_with_COMPASS_Subsets/"))

# Arial font setup. Downloaded afms from https://github.com/microsoft/microsoft-r-open/tree/ec3fd89e5fb5794bd8149905c134ad801bb61800
Arial <- Type1Font(family = "Arial",
                   metrics = c(here::here("data/Arial_afm/ArialMT.afm"),
                               here::here("data/Arial_afm/ArialMT-Bold.afm"),
                               here::here("data/Arial_afm/ArialMT-Italic.afm"),
                               here::here("data/Arial_afm/ArialMT-BoldItalic.afm")))
windowsFonts("Arial" = windowsFont("Arial"))
pdfFonts(Arial = Arial)

# Extract frequencies
nodes <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/Naive",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/TCM",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/TEM",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/TEMRA",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/IFNg+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/IL17a+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/TBET+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/RORyT+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/CXCR3+CCR6+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/CXCR3+CCR6-",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4 Positive/CD4+/CD4_COMPASS_Subsets/CXCR3-CCR6+")
nodes_short <- str_replace(nodes, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/CD4 Positive\\/CD4\\+\\/CD4_COMPASS_Subsets\\/", "")

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

# Add column indicating "Node" in memory and Th marker count dataframes
naive_counts <- count_list[[nodes_short[1]]]
naive_counts["Node"] <- nodes_short[1]

tcm_counts <- count_list[[nodes_short[2]]] 
tcm_counts["Node"] <- nodes_short[2]

tem_counts <- count_list[[nodes_short[3]]] 
tem_counts["Node"] <- nodes_short[3]

temra_counts <- count_list[[nodes_short[4]]] 
temra_counts["Node"] <- nodes_short[4]

cxcr3_ccr6_counts <- count_list[[nodes_short[9]]] 
cxcr3_ccr6_counts["Node"] <- nodes_short[9]

cxcr3_counts <- count_list[[nodes_short[10]]] 
cxcr3_counts["Node"] <- nodes_short[10]

ccr6_counts <- count_list[[nodes_short[11]]] 
ccr6_counts["Node"] <- nodes_short[11]

# Combine memory count dataframes and Th marker count dataframes
mem_counts <- bind_rows(naive_counts, tcm_counts) %>%
  bind_rows(tem_counts) %>%
  bind_rows(temra_counts)

mem_nodes_compare <- c(nodes_short[1], nodes_short[2], nodes_short[3], nodes_short[4])

th_counts <- bind_rows(cxcr3_ccr6_counts, cxcr3_counts) %>%
  bind_rows(ccr6_counts) 

th_nodes_compare <- c(nodes_short[9], nodes_short[10], nodes_short[11])

# Set color scheme and stims
fill_colors <- c("Pneg" = "#984EA3", "TST+" = "#4DAF4A")
stims <- c("DMSO", "PP1", "TB WCL")

# Get plots
naive_plots <- purrr::pmap(.l = list(stims),
                         .f = function(n) {
                           make_mag_plots(count_list[[nodes_short[1]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                          paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[1],
                                          y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 60))
                         })
names(naive_plots) <- stims

tcm_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[2]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[2],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(30, 100))
                              })
names(tcm_plots) <- stims

tem_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[3]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[3],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 40))
                                })
names(tem_plots) <- stims

temra_plots <- purrr::pmap(.l = list(stims),
                         .f = function(n) {
                           make_mag_plots(count_list[[nodes_short[4]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                          paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[4],
                                          y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 1.5))
                         })
names(temra_plots) <- stims

mem_plots <- purrr::pmap(.l = list(rep(stims, each = 2), rep(c("Pneg", "TST+"), times = 3)),
                          .f = function(n, status) {
                            multi_mag_plot(counts = mem_counts, stim_keep = n, status_keep = status, fill_colors = fill_colors, nodes_compare = mem_nodes_compare, ylim = c(0, 90))
                          })
names(mem_plots) <- paste0(rep(stims, each = 2), "_", rep(c("Pneg", "TST+"), times = 3))

ifng_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[5]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[5],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 10))
                                })
names(ifng_plots) <- stims

il17a_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[6]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[6],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 0.6))
                                })
names(il17a_plots) <- stims

tbet_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[7]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[7],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 75))
                                })
names(tbet_plots) <- stims

roryt_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[8]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[8],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 25))
                                })
names(roryt_plots) <- stims

cxcr3_ccr6_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[9]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[9],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 40))
                                })
names(cxcr3_ccr6_plots) <- stims

cxcr3_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[10]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[10],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 40))
                                })
names(cxcr3_plots) <- stims

ccr6_plots <- purrr::pmap(.l = list(stims),
                           .f = function(n) {
                             make_mag_plots(count_list[[nodes_short[11]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                            paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[11],
                                            y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 40))
                           })
names(ccr6_plots) <- stims

th_plots <- purrr::pmap(.l = list(rep(stims, each = 2), rep(c("Pneg", "TST+"), times = 3)),
                          .f = function(n, status) {
                            multi_mag_plot(counts = th_counts, stim_keep = n, status_keep = status, fill_colors = fill_colors, nodes_compare = th_nodes_compare, ylim = c(0, 45))
                          })
names(th_plots) <- paste0(rep(stims, each = 2), "_", rep(c("Pneg", "TST+"), times = 3))

# Save plots
cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_Naive_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(naive_plots[["DMSO"]])
print(naive_plots[["PP1"]])
print(naive_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_TCM_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(tcm_plots[["DMSO"]])
print(tcm_plots[["PP1"]])
print(tcm_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_TEM_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(tem_plots[["DMSO"]])
print(tem_plots[["PP1"]])
print(tem_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_TEMRA_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(temra_plots[["DMSO"]])
print(temra_plots[["PP1"]])
print(temra_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_Mem_multi_plots.pdf"), width=6, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(mem_plots[["DMSO_Pneg"]])
print(mem_plots[["DMSO_TST+"]])
print(mem_plots[["PP1_Pneg"]])
print(mem_plots[["PP1_TST+"]])
print(mem_plots[["TB WCL_Pneg"]])
print(mem_plots[["TB WCL_TST+"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_IFNg_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(ifng_plots[["DMSO"]])
print(ifng_plots[["PP1"]])
print(ifng_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_IL17a_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(il17a_plots[["DMSO"]])
print(il17a_plots[["PP1"]])
print(il17a_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_TBET_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(tbet_plots[["DMSO"]])
print(tbet_plots[["PP1"]])
print(tbet_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_RORyT_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(roryt_plots[["DMSO"]])
print(roryt_plots[["PP1"]])
print(roryt_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_CXCR3_CCR6_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cxcr3_ccr6_plots[["DMSO"]])
print(cxcr3_ccr6_plots[["PP1"]])
print(cxcr3_ccr6_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_CXCR3_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cxcr3_plots[["DMSO"]])
print(cxcr3_plots[["PP1"]])
print(cxcr3_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_CCR6_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(ccr6_plots[["DMSO"]])
print(ccr6_plots[["PP1"]])
print(ccr6_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Th_CD4_COMPASS_Subsets_Analysis/Th_CD4_COMPASS_Th_multi_plots.pdf"), width=6, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(th_plots[["DMSO_Pneg"]])
print(th_plots[["DMSO_TST+"]])
print(th_plots[["PP1_Pneg"]])
print(th_plots[["PP1_TST+"]])
print(th_plots[["TB WCL_Pneg"]])
print(th_plots[["TB WCL_TST+"]])
dev.off()
