# Individual magnitude/frequency plots

make_mag_plots <- function(counts, counts_no_outlier = NULL, current_stim,
                           num_comparisons, groups_to_compare, paired, adjust_p, fill_colors,
                           group_by_colname, subtitle, ylim = NULL, y_axis_text,
                           y_axis_size, facet_by = NULL) {
  
  counts <- counts %>%
    dplyr::filter(Stim == current_stim)
  
  if(!is.null(counts_no_outlier)) {
    counts_to_plot <- counts_no_outlier %>%
        dplyr::filter(Stim == current_stim)
  } else {
    counts_to_plot <- counts
  }
  
  # P-values are unadjusted
  fmla <- formula(paste("Freq ~ ", group_by_colname))
  
  if(paired) {
    # Signed-rank test
    test <- wilcox.test(fmla, data = counts, paired = TRUE)
  } else {
    # Rank sum test
    test <- wilcox.test(fmla, data = counts, paired = FALSE)
  }
 
  if(adjust_p) {
    test_df <- data.frame(p = as.numeric(unlist(test)["p.value"])) %>%
      mutate(p.adj = p.adjust(p, method = "bonferroni", n = num_comparisons)) %>% # Strict
      mutate(p_val_text = if_else(p.adj < 0.001, "p<0.001", paste0("p=", formatC(round(p.adj, 3), format='f', digits=3))))
  } else {
    test_df <- data.frame(p = as.numeric(unlist(test)["p.value"])) %>%
      mutate(p_val_text = if_else(p < 0.001, "p<0.001", paste0("p=", formatC(round(p, 3), format='f', digits=3))))
  }
  
  current_plot <- ggplot(counts_to_plot, aes(x = !!as.symbol(group_by_colname), y = Freq)) +
    geom_boxplot(outlier.shape=NA, position = position_dodge2(preserve = "total")) +
    geom_quasirandom(size=3, shape = 16, width = 0.3, aes(color=!!as.symbol(group_by_colname))) +
    scale_color_manual(values = fill_colors) +
    scale_x_discrete(labels = stringr::str_replace_all(groups_to_compare, " ", "\n"), expand = c(0.3, 0.3))
  
  current_plot <- current_plot +
    theme_bw() +
    theme(text = element_text(family="Arial"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = y_axis_size),
          axis.text.y = element_text(color="black", size=17),
          axis.text.x = element_text(color="black", size=20),
          plot.title = element_text(hjust = 0.5, size=21),
          plot.subtitle = element_text(hjust = 0.5, size=13),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(0.3, 0.2, 0.1, 0.2, "cm")) +
    labs(title = as.character(current_stim),
         subtitle = subtitle,
         y = y_axis_text) +
    force_panelsizes(rows = unit(3.5, "in"),
                     cols = unit(3, "in"))
  
  if(!is.null(facet_by)) {
    fmla2 <- formula(paste(". ~ ", facet_by))
    current_plot <- current_plot +
      facet_grid(fmla2) +
      stat_compare_means(comparisons = list(groups_to_compare), label = "p.format",
                         method = "wilcox.test", paired = FALSE, tip.length = 0)
  } else {
    if(!is.null(ylim)) {
      current_plot <- current_plot + 
        coord_cartesian(ylim = ylim) +
        annotate("text", x = 1.5, y = ylim[2] - 0.09*diff(ylim),
                 label = test_df$p_val_text, size=5.5)
    } else {
      plot_ylims <- ggplot_build(current_plot)$layout$panel_params[[1]]$y.range
      current_plot <- current_plot + 
        annotate("text", x = 1.5, y = plot_ylims[2] + 0.01*diff(plot_ylims),
                 label = test_df$p_val_text, size=5.5) +
        coord_cartesian(ylim = c(plot_ylims[[1]], plot_ylims[[2]] + 0.09*diff(plot_ylims)))
    } 
  }

}

############################################################################################################################

multi_mag_plot <- function(counts, current_status = "Pneg", fill_colors = NULL, subtitle, groups_to_compare, y_axis_text) {
  counts <- counts %>% 
    dplyr::filter(Status == current_status)
  
  counts_wide <- counts %>%
    pivot_wider(id_cols = c(`SAMPLE ID`, Status), names_from = Stim, values_from = Freq) %>%
    na.omit() # can't perform signed-rank test with missing data

  # Signed-rank test
  one_two_signed_rank_result <- wilcox.test(counts_wide$DMSO, counts_wide$PP1, paired = T)
  two_three_signed_rank_result <- wilcox.test(counts_wide$PP1, counts_wide$`TB WCL`, paired = T)
  one_three_signed_rank_result <- wilcox.test(counts_wide$DMSO, counts_wide$`TB WCL`, paired = T)

  stim_x_order <- 1:3
  counts$Stim <- factor(counts$Stim, levels = groups_to_compare)
  names(stim_x_order) <- levels(counts$Stim)
  test_df <- data.frame(group1 = c("DMSO", "PP1", "DMSO"),
                        group2 = c("PP1", "TB WCL", "TB WCL"),
                        p.val = c(one_two_signed_rank_result$p.value,
                                  two_three_signed_rank_result$p.value,
                                  one_three_signed_rank_result$p.value)) %>%
    # mutate(p.adj = p.adjust(p.val, method = "bonferroni")) %>% # Strict
    mutate(group1.xloc = unname(stim_x_order[group1]),
           group2.xloc = unname(stim_x_order[group2]),
           p_val_text = sapply(p.val, function(p) {
             if(p < 0.001) {
               "p<0.001"
             } else {
               paste0("p=", round(p, 3))
             }
           }),
           y_pos = c(counts %>% dplyr::filter(Stim %in% c("DMSO", "PP1")) %>%
                       dplyr::pull(Freq) %>% max() + 0.1*diff(range(counts %>% dplyr::pull(Freq))),
                     counts %>% dplyr::filter(Stim %in% c("PP1", "TB WCL")) %>%
                       dplyr::pull(Freq) %>% max() + 0.1*diff(range(counts %>% dplyr::pull(Freq))),
                     counts %>% dplyr::filter(Stim %in% c("DMSO", "PP1", "TB WCL")) %>%
                       dplyr::pull(Freq) %>% max() + 0.25*diff(range(counts %>% dplyr::pull(Freq)))),
           geom_signif_group = paste0(group1, group2))

  medians <- counts %>% group_by(Stim) %>%
    summarise(med = median(Freq)) %>% ungroup() %>%
    mutate(x.min.segment = 1:3 - 0.1,
           x.end.segment = 1:3 + 0.1)
  
  current_plot <- ggplot(counts, aes(Stim, Freq)) +
    geom_boxplot(outlier.shape=NA, position = position_dodge2(preserve = "total")) +
    geom_quasirandom(size=3, shape = 16, width = 0.3, aes(color=!!as.symbol("Status"))) +
    labs(y = y_axis_text,
         title = current_status,
         subtitle = subtitle) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(color="black", size=17),
          axis.text.x = element_text(color="black", size=17),
          plot.title = element_text(hjust = 0.5, size=21),
          plot.subtitle = element_text(hjust = 0.5, size=13),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(0.3, 0.2, 0.1, 0.2, "cm")) +
    scale_x_discrete(expand = c(0.15, 0.15),
                     labels = groups_to_compare) +
    scale_color_manual(values = fill_colors[[current_status]]) +
    ggsignif::geom_signif(inherit.aes=F,data=test_df,
                          aes_string(xmin="group1.xloc", xmax="group2.xloc",
                                     annotations="p_val_text", y_position="y_pos",
                                     group="geom_signif_group"), 
                          tip_length = c(0, 0),
                          textsize=5,
                          size = 0.75,
                          manual = TRUE) +
    coord_cartesian(ylim = c(NA, max(test_df$y_pos) + 0.1*diff(range(counts %>% dplyr::pull(Freq)))))
}
