# This script contains:
# 1) All functions in about.Rmd
# 2) Any global variables (such as the dataframe of LowN response color codes)
# to be sourced in other scripts

#### Environment ####
sapply(c("dplyr", "readr", "tidyr", "purrr", "ggplot2", "ggpubr", "matrixStats"), require, character.only=TRUE)
load("data_files/CleanedCounts.RData")

#### Global Variables ####
# colordf
response_colordf <- tibble(limits = c("none", "early", "late", "gradual", "both"),
                           value = c("#A1A3A6", "#000000", "#000000", "#000000", "#000000"))
conservation_colordf <- tibble(limits = c("conserved", "Scer", "Spar", "Hyb", "ScerHyb", "SparHyb", "ScerSpar", "unresponsive"),
                           value = c("#FFD540", "#F15F5B", "#47AA48", "#28A5DE", "#7D70B3", "#18B682", "#E4E5E5", "#AAACAF"))
# shorthand b/c we use these a lot
SparColor <- conservation_colordf$value[conservation_colordf$limits == "Spar"]
ScerColor <- conservation_colordf$value[conservation_colordf$limits == "Scer"]

# goslim
# table of gene ontology terms associated with each gene
# each row is a term, so this is a very long table
goslim <- read.table("data_files/downloaded_genomes_and_features/go_slim_mapping.tab", header = FALSE, sep = "\t") |>
  as_tibble()
colnames(goslim) <- c("ORF", # (mandatory) - Systematic name of the gene (or gene complex if it starts with CPX-)
                      "gene", # (optional) - Gene name, if one exists
                      "SGDID", # SGDID (mandatory) - the SGDID, unique database identifier for the gene
                      "GO_aspect", # (mandatory) - which ontology: P=Process, F=Function, C=Component
                      "GOslim_term", # (mandatory) - the name of the GO term that was selected as a GO Slim term
                      "GOID", # (optional) - the unique numerical identifier of the GO term
                      "feature_type") # (mandatory) - a description of the sequence feature, such as ORF or tRNA
too_vague_terms <- c("cellular process", "molecular function", "biological process")
goslim <- filter(goslim, !(GOslim_term %in% too_vague_terms) & 
                   (ORF %in% rownames(counts)))

#### Functions ####
# plots WT vs TFdel in both species (option to exclude WT or TFdel)
plotExpressionProfileTFdel <- function(.cts1, .cts2, .cts3, .cts4,
                                       .info1, .info2, .info3, .info4,
                                       .name1 = "S. cerevisiae WT",
                                       .name2 = "S. paradoxus WT",
                                       .name3 = "S. cerevisiae TFdel",
                                       .name4 = "S. paradoxus TFdel",
                                       .color1 = ScerColor,
                                       .color2 = SparColor,
                                       .color3 = ScerColor,
                                       .color4 = SparColor,
                                       .normalization = c("none", "log2", "scale", "centered log2"),
                                       .show_points_wt = TRUE,
                                       .show_points_tfdel = TRUE,
                                       .show_lines_wt = TRUE,
                                       .show_lines_tfdel = FALSE,
                                       .show_wt = TRUE,
                                       .show_tfdel = TRUE,
                                       .show_confidence_intervals = TRUE,
                                       .plotlims = NULL) {
  if (.normalization == "none") {
    norm_func <- identity
    ylabel <- "Expression"
  }
  if (.normalization == "log2") {
    norm_func <- \(x) {log2(x + 1)}
    ylabel <- "Expression (log2)"
  }
  if (.normalization == "scale") {
    norm_func <- \(x) {t(scale(t(x)))}
    ylabel <- "Expression (centered and scaled)"
  }
  if (.normalization == "centered log2") {
    norm_func <- \(x) {log2(x + 1) - rowMeans(log2(x + 1), na.rm = TRUE)}
    ylabel <- "Expression (centered log2)"
  }
  if (!setequal(unique(.info1$experiment), unique(.info2$experiment))) {
    stop("sample info dataframes do not contain same set of experiments\n")
  }
  info1 <- tibble(experiment = "LowN",
                  time_point_str = .info1$time_point_str,
                  genotype = .info1$genotype,
                  well_flask_ID = .info1$well_flask_ID)
  info2 <- tibble(experiment = "LowN",
                  time_point_str = .info2$time_point_str,
                  genotype = .info2$genotype,
                  well_flask_ID = .info2$well_flask_ID)
  info3 <- tibble(experiment = "LowN",
                  time_point_str = .info3$time_point_str,
                  genotype = .info3$genotype,
                  well_flask_ID = .info3$well_flask_ID)
  info4 <- tibble(experiment = "LowN",
                  time_point_str = .info4$time_point_str,
                  genotype = .info4$genotype,
                  well_flask_ID = .info4$well_flask_ID)
  expr1 <- norm_func(.cts1) |> t()
  expr2 <- norm_func(.cts2) |> t()
  expr3 <- norm_func(.cts3) |> t()
  expr4 <- norm_func(.cts4) |> t()
  colnames(expr1) <- rownames(.cts1)
  colnames(expr2) <- rownames(.cts2)
  colnames(expr3) <- rownames(.cts3)
  colnames(expr4) <- rownames(.cts4)
  gdf1 <- bind_cols(expr1, info1) |> 
    pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
  gdf1$group_id <- "1"
  gdf2 <- bind_cols(expr2, info2) |> 
    pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
  gdf2$group_id <- "2"
  gdf3 <- bind_cols(expr3, info3) |> 
    pivot_longer(cols = colnames(expr3), names_to = "gene_name", values_to = "expr")
  gdf3$group_id <- "3"
  gdf4 <- bind_cols(expr4, info4) |> 
    pivot_longer(cols = colnames(expr4), names_to = "gene_name", values_to = "expr")
  gdf4$group_id <- "4"
  # converting each gene's expression to its mean expression between replicates
  gdf <- bind_rows(gdf1, gdf2, gdf3, gdf4) |> 
    group_by(group_id, time_point_str, well_flask_ID, genotype) |> 
    summarise(expr = mean(expr, na.rm = TRUE)) |> ungroup()
  plotlimdf <- gdf |> group_by(time_point_str, genotype, group_id) |>
    summarise(mean_expr = mean(expr, na.rm = TRUE),
              sd_expr = sd(expr, na.rm = TRUE)) 
  max_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> max(na.rm = TRUE)
  min_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> min(na.rm = TRUE)
  buffer <- 0.25
  max_avg <- max_avg + buffer
  min_avg <- min_avg - buffer
  # min/maxs when plotting points as well as averages:
  max_expr <- max(gdf$expr, na.rm = TRUE)
  min_expr <- min(gdf$expr, na.rm = TRUE)
  # setting plotlims if they haven't been manually set
  if (is.null(.plotlims)) {
    if (!(.show_points_wt | .show_points_tfdel)) {
      .plotlims <- c(min_avg, max_avg)
    }
    if (.show_points_wt | .show_points_tfdel) {
      .plotlims <- c(min_expr - 0.1, max_expr + 0.1) # buffer b/c TFdel points are very large
      if (.normalization == "none" & min_expr == 0) { # for some reason, points right at 0 get cut off, probably the jitter is too significant without normalization
        .plotlims <- c(min_expr - 5, max_expr + 0.1)
      }
    }
  }
  plotdf <- gdf
  # plotting
  p <- ggplot() + 
    theme_classic() +
    scale_color_discrete(type = c(.color1, .color2, .color3, .color4), 
                         labels = c(.name1, .name2, .name3, .name4),
                         limits = c("1", "2", "3", "4")) +
    theme(legend.title = element_blank(),
          legend.position = "none") +
    ylab("") +
    xlab("") +
    ylim(.plotlims)
  if (.show_points_wt & .show_wt) {
    p <- p + geom_jitter(data = filter(plotdf, genotype == "WT"),
                         aes(x = time_point_str, y = expr, color = group_id), 
                         size = 0.3, alpha = 1, shape = ".") 
  }
  if (.show_lines_wt & .show_wt) {
    # lines trace average expr at each timepoint/experiment for each group
    avgexpr1 <- plotdf %>% filter(group_id == 1) |> group_by(time_point_str, group_id) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                    sd_expr = sd(expr, na.rm = TRUE))
    avgexpr2 <- plotdf %>% filter(group_id == 2) |> group_by(time_point_str, group_id) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                    sd_expr = sd(expr, na.rm = TRUE))
    # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
    p <- p +
      geom_line(data = avgexpr1, aes(x = time_point_str, y = mean_expr, group = group_id), color = .color1, linewidth = 1, linetype = "solid") +
      geom_line(data = avgexpr2, aes(x = time_point_str, y = mean_expr, group = group_id), color = .color2, linewidth = 1, linetype = "solid")
    if (.show_confidence_intervals) {
      # calculating 95% confidence in the mean
      nGenes <- nrow(.cts1)
      avgexpr1$CI <- 1.96*(avgexpr1$sd_expr/sqrt(nGenes))
      avgexpr2$CI <- 1.96*(avgexpr2$sd_expr/sqrt(nGenes))
      p <- p + 
        geom_ribbon(data = avgexpr1, aes(x = time_point_str, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                    fill = .color1, alpha = 0.3) +
        geom_ribbon(data = avgexpr2, aes(x = time_point_str, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                    fill = .color2, alpha = 0.3)
    }
    if (.show_points_tfdel & .show_tfdel) {
      p <- p + geom_jitter(data = filter(plotdf, genotype != "WT"),
                           aes(x = time_point_str, y = expr, color = group_id),
                           shape = "+", size = 6, width = 0.05, height = 0) # height = 0 prevents vertical jitter
    }
    if (.show_lines_tfdel & .show_tfdel) {
      # lines trace average expr at each timepoint/experiment for each group
      avgexpr3 <- plotdf %>% filter(group_id == 3) |> group_by(time_point_str, group_id) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                      sd_expr = sd(expr, na.rm = TRUE))
      avgexpr4 <- plotdf %>% filter(group_id == 4) |> group_by(time_point_str, group_id) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                      sd_expr = sd(expr, na.rm = TRUE))
      # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
      p <- p +
        geom_line(data = avgexpr3, aes(x = time_point_str, y = mean_expr, group = group_id), color = .color3, linewidth = 1, linetype = "dashed") +
        geom_line(data = avgexpr4, aes(x = time_point_str, y = mean_expr, group = group_id), color = .color4, linewidth = 1, linetype = "dashed")
      if (.show_confidence_intervals) {
        # calculating 95% confidence in the mean
        nGenes <- nrow(.cts1)
        avgexpr3$CI <- 1.96*(avgexpr3$sd_expr/sqrt(nGenes))
        avgexpr4$CI <- 1.96*(avgexpr4$sd_expr/sqrt(nGenes))
        p <- p + 
          geom_ribbon(data = avgexpr3, aes(x = time_point_str, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                      fill = .color1, alpha = 0.3) +
          geom_ribbon(data = avgexpr4, aes(x = time_point_str, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                      fill = .color2, alpha = 0.3)
      }
    }
    return(p)
  }
}

# Plots each genotype as a separate line
plotExpressionProfileTFdels <- function(.cts,
                                        .info,
                                        .name1 = "WT",
                                        .name2 = "TFdel",
                                        .color1 = "red",
                                        .color2 = "grey",
                                        .normalization = c("none", "log2", "scale", "centered log2"),
                                        .ylims = NULL,
                                        .ylab = TRUE) {
  if (.normalization == "none") {
    norm_func <- identity
    ylabel <- "Expression"
  }
  if (.normalization == "log2") {
    norm_func <- \(x) {log2(x + 1)}
    ylabel <- "Expression (log2)"
  }
  if (.normalization == "scale") {
    norm_func <- scale
    ylabel <- "Expression (centered and scaled)"
  }
  if (.normalization == "centered log2") {
    norm_func <- \(x) {log2(x + 1) - rowMeans(log2(x + 1), na.rm = TRUE)}
    ylabel <- "Expression (centered log2)"
  }
  if (!.ylab) {
    ylabel <- ""
  }
  gdf <- bind_cols(tibble(sample_name = colnames(.cts)),
                   t(.cts)) |> 
    left_join(.info, by = "sample_name") |> 
    pivot_longer(cols = rownames(.cts), 
                 names_to = "gene_name", values_to = "expr")
  plotdf <- gdf |> group_by(genotype, gene_name) |> 
    reframe(expr = norm_func(expr),
            time_point_str = time_point_str) |>
    ungroup() |>
    group_by(genotype, time_point_str) |> 
    summarise(mean_expr = mean(expr, na.rm = TRUE))
  # plotting
  p <- ggplot(plotdf, aes(x = time_point_str, y = mean_expr)) +
    geom_line(data = filter(plotdf, genotype != "WT"),
              aes(group = genotype), color = .color2, alpha = 0.5) +
    geom_line(data = filter(plotdf, genotype == "WT"),
              aes(group = genotype), color = .color1,
              linewidth = 1.5) +
    theme_classic() +
    ylab(ylabel) +
    xlab("")
  if (!is.null(.ylims)) {
    p <- p + ylim(.ylims)
  }
  return(p)
}

# This function is in need of description, but my guess is that it treats every gene as its own line rather than plotting the mean of all genes in the group
plotSingleProfilesTFdel <- function(.cts1, .cts2, .cts3, .cts4,
                                    .info1, .info2, .info3, .info4,
                                    .name1 = "S. cerevisiae WT",
                                    .name2 = "S. paradoxus WT",
                                    .name3 = "S. cerevisiae TFdel",
                                    .name4 = "S. paradoxus TFdel",
                                    .color1 = "orange1",
                                    .color2 = "blue2",
                                    .color3 = "orange4",
                                    .color4 = "blue4",
                                    .normalization = c("none", "log2", "scale", "centered log2"),
                                    .show_points_wt = TRUE,
                                    .show_points_tfdel = TRUE,
                                    .show_lines_wt = TRUE,
                                    .show_lines_tfdel = FALSE,
                                    .show_wt = TRUE,
                                    .show_tfdel = TRUE,
                                    .show_confidence_intervals = TRUE,
                                    .plotlims = NULL) {
  if (.normalization == "none") {
    norm_func <- identity
    ylabel <- "Expression"
  }
  if (.normalization == "log2") {
    norm_func <- \(x) {log2(x + 1)}
    ylabel <- "Expression (log2)"
  }
  if (.normalization == "scale") {
    norm_func <- \(x) {t(scale(t(x)))}
    ylabel <- "Expression (centered and scaled)"
  }
  if (.normalization == "centered log2") {
    norm_func <- \(x) {log2(x + 1) - rowMeans(log2(x + 1), na.rm = TRUE)}
    ylabel <- "Expression (centered log2)"
  }
  if (!setequal(unique(.info1$experiment), unique(.info2$experiment))) {
    stop("sample info dataframes do not contain same set of experiments\n")
  }
  info1 <- tibble(experiment = "LowN",
                  time_point_str = .info1$time_point_str,
                  genotype = .info1$genotype,
                  well_flask_ID = .info1$well_flask_ID)
  info2 <- tibble(experiment = "LowN",
                  time_point_str = .info2$time_point_str,
                  genotype = .info2$genotype,
                  well_flask_ID = .info2$well_flask_ID)
  info3 <- tibble(experiment = "LowN",
                  time_point_str = .info3$time_point_str,
                  genotype = .info3$genotype,
                  well_flask_ID = .info3$well_flask_ID)
  info4 <- tibble(experiment = "LowN",
                  time_point_str = .info4$time_point_str,
                  genotype = .info4$genotype,
                  well_flask_ID = .info4$well_flask_ID)
  expr1 <- norm_func(.cts1) |> t()
  expr2 <- norm_func(.cts2) |> t()
  expr3 <- norm_func(.cts3) |> t()
  expr4 <- norm_func(.cts4) |> t()
  colnames(expr1) <- rownames(.cts1)
  colnames(expr2) <- rownames(.cts2)
  colnames(expr3) <- rownames(.cts3)
  colnames(expr4) <- rownames(.cts4)
  gdf1 <- bind_cols(expr1, info1) |> 
    pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
  gdf1$group_id <- "1"
  gdf2 <- bind_cols(expr2, info2) |> 
    pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
  gdf2$group_id <- "2"
  gdf3 <- bind_cols(expr3, info3) |> 
    pivot_longer(cols = colnames(expr3), names_to = "gene_name", values_to = "expr")
  gdf3$group_id <- "3"
  gdf4 <- bind_cols(expr4, info4) |> 
    pivot_longer(cols = colnames(expr4), names_to = "gene_name", values_to = "expr")
  gdf4$group_id <- "4"
  # enumerating replicates, so each one can stack right on top of each other
  plotdf <- bind_rows(gdf1, gdf2, gdf3, gdf4)
  plotdf <- plotdf |> group_by(group_id, genotype, time_point_str, gene_name) |> 
    reframe(rep = rank(gene_name, ties.method = "first"),
            gene_name = gene_name,
            time_point_str = time_point_str,
            genotype = genotype,
            expr = expr,
            group_id = group_id)
  plotdf$rep <- if_else(plotdf$genotype == "WT",
                        true = 1, false = plotdf$rep)
  # plotting
  p <- ggplot() + 
    theme_classic() +
    scale_color_discrete(type = c(.color1, .color2, .color3, .color4), 
                         labels = c(.name1, .name2, .name3, .name4),
                         limits = c("1", "2", "3", "4")) +
    theme(legend.title = element_blank(),
          legend.position = "none") +
    ylab("") +
    xlab("")
  # lines trace average expr at each timepoint/experiment for each group
  avgexpr1 <- plotdf %>% filter(group_id == 1) |> group_by(time_point_str, group_id, gene_name) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                             sd_expr = sd(expr, na.rm = TRUE))
  avgexpr2 <- plotdf %>% filter(group_id == 2) |> group_by(time_point_str, group_id, gene_name) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                             sd_expr = sd(expr, na.rm = TRUE))
  p <- p +
    geom_line(data = avgexpr1, aes(x = time_point_str, y = mean_expr, group = interaction(group_id, gene_name)), color = .color1) +
    geom_line(data = avgexpr2, aes(x = time_point_str, y = mean_expr, group = interaction(group_id, gene_name)), color = .color2) +
    geom_line(data =  filter(plotdf, group_id == 3),
              position = position_nudge(x = 0.05),
              aes(x = time_point_str, y = expr,
                  color = group_id, group = interaction(group_id, time_point_str, gene_name))) +
    geom_point(data =  filter(plotdf, group_id == 3),
               position = position_nudge(x = 0.05),
               size = 0.25,
               aes(x = time_point_str, y = expr,
                   color = group_id)) +
    geom_line(data =  filter(plotdf, group_id == 4),
              position = position_nudge(x = -0.05),
              aes(x = time_point_str, y = expr,
                  color = group_id, group = interaction(group_id, time_point_str, gene_name))) +
    geom_point(data =  filter(plotdf, group_id == 4),
               position = position_nudge(x = -0.05),
               size = 0.25,
               aes(x = time_point_str, y = expr,
                   color = group_id)) +
    facet_wrap(~gene_name)
  return(p)
}

# Wrapper function for plotExpressionProfileTFdel and plotSingleProfilesTFdel
plotGenesTFdel <- function(.gene_idxs, .tf, .parents_or_hybrid = "parents",
                           .plotlims = NULL, .single_genes = FALSE,
                           .normalization = "log2", .show_wt = TRUE, .show_tfdel = TRUE) {
  if (.parents_or_hybrid == "parents" & !.single_genes) {
    p <- plotExpressionProfileTFdel(.cts1 = counts[.gene_idxs,
                                                   sample_info$organism == "cer" &
                                                     sample_info$genotype == "WT",
                                                   drop = FALSE],
                                    .cts2 = counts[.gene_idxs,
                                                   sample_info$organism == "par" &
                                                     sample_info$genotype == "WT",
                                                   drop = FALSE],
                                    .cts3 = counts[.gene_idxs,
                                                   sample_info$organism == "cer" &
                                                     sample_info$genotype == paste0(.tf, "delete"),
                                                   drop = FALSE],
                                    .cts4 = counts[.gene_idxs,
                                                   sample_info$organism == "par" &
                                                     sample_info$genotype == paste0(.tf, "delete"),
                                                   drop = FALSE],
                                    .info1 = sample_info[sample_info$organism == "cer" &
                                                           sample_info$genotype == "WT",],
                                    .info2 = sample_info[sample_info$organism == "par" &
                                                           sample_info$genotype == "WT",],
                                    .info3 = sample_info[sample_info$organism == "cer" &
                                                           sample_info$genotype == paste0(.tf, "delete"),],
                                    .info4 = sample_info[sample_info$organism == "par" &
                                                           sample_info$genotype == paste0(.tf, "delete"),],
                                    .normalization = .normalization, .plotlims = .plotlims,
                                    .show_wt = .show_wt, .show_tfdel = .show_tfdel)
  }
  if (.parents_or_hybrid == "hybrid" & !.single_genes) {
    p <- plotExpressionProfileTFdel(.cts1 = counts[.gene_idxs,
                                                   sample_info$organism == "hyb" &
                                                     sample_info$allele == "cer" &
                                                     sample_info$genotype == "WT",
                                                   drop = FALSE],
                                    .cts2 = counts[.gene_idxs, 
                                                   sample_info$organism == "hyb" &
                                                     sample_info$allele == "par" &
                                                     sample_info$genotype == "WT",
                                                   drop = FALSE],
                                    .cts3 = counts[.gene_idxs,
                                                   sample_info$organism == "hyb" &
                                                     sample_info$allele == "cer" &
                                                     sample_info$genotype == paste0(.tf, "delete"),
                                                   drop = FALSE],
                                    .cts4 = counts[.gene_idxs,
                                                   sample_info$organism == "hyb" &
                                                     sample_info$allele == "par" &
                                                     sample_info$genotype == paste0(.tf, "delete"),
                                                   drop = FALSE],
                                    .info1 = sample_info[sample_info$organism == "hyb" &
                                                           sample_info$allele == "cer" &
                                                           sample_info$genotype == "WT",],
                                    .info2 = sample_info[sample_info$organism == "hyb" &
                                                           sample_info$allele == "par" &
                                                           sample_info$genotype == "WT",],
                                    .info3 = sample_info[sample_info$organism == "hyb" &
                                                           sample_info$allele == "cer" &
                                                           sample_info$genotype == paste0(.tf, "delete"),],
                                    .info4 = sample_info[sample_info$organism == "hyb" &
                                                           sample_info$allele == "par" &
                                                           sample_info$genotype == paste0(.tf, "delete"),],
                                    .normalization = .normalization, .plotlims = .plotlims,
                                    .show_wt = .show_wt, .show_tfdel = .show_tfdel)
  }
  if (.parents_or_hybrid == "parents" & .single_genes) {
    p <- plotSingleProfilesTFdel(.cts1 = counts[.gene_idxs,
                                                sample_info$organism == "cer" &
                                                  sample_info$genotype == "WT",
                                                drop = FALSE],
                                 .cts2 = counts[.gene_idxs,
                                                sample_info$organism == "par" &
                                                  sample_info$genotype == "WT",
                                                drop = FALSE],
                                 .cts3 = counts[.gene_idxs,
                                                sample_info$organism == "cer" &
                                                  sample_info$genotype == paste0(.tf, "delete"),
                                                drop = FALSE],
                                 .cts4 = counts[.gene_idxs,
                                                sample_info$organism == "par" &
                                                  sample_info$genotype == paste0(.tf, "delete"),
                                                drop = FALSE],
                                 .info1 = sample_info[sample_info$organism == "cer" &
                                                        sample_info$genotype == "WT",],
                                 .info2 = sample_info[sample_info$organism == "par" &
                                                        sample_info$genotype == "WT",],
                                 .info3 = sample_info[sample_info$organism == "cer" &
                                                        sample_info$genotype == paste0(.tf, "delete"),],
                                 .info4 = sample_info[sample_info$organism == "par" &
                                                        sample_info$genotype == paste0(.tf, "delete"),],
                                 .normalization = .normalization, .plotlims = .plotlims,
                                 .show_wt = .show_wt, .show_tfdel = .show_tfdel)
  }
  if (.parents_or_hybrid == "hybrid" & .single_genes) {
    p <- plotSingleProfilesTFdel(.cts1 = counts[.gene_idxs,
                                                sample_info$organism == "hyb" &
                                                  sample_info$allele == "cer" &
                                                  sample_info$genotype == "WT",
                                                drop = FALSE],
                                 .cts2 = counts[.gene_idxs,
                                                sample_info$organism == "hyb" &
                                                  sample_info$allele == "par" &
                                                  sample_info$genotype == "WT",
                                                drop = FALSE],
                                 .cts3 = counts[.gene_idxs,
                                                sample_info$organism == "hyb" &
                                                  sample_info$allele == "cer" &
                                                  sample_info$genotype == paste0(.tf, "delete"),
                                                drop = FALSE],
                                 .cts4 = counts[.gene_idxs,
                                                sample_info$organism == "hyb" &
                                                  sample_info$allele == "par" &
                                                  sample_info$genotype == paste0(.tf, "delete"),
                                                drop = FALSE],
                                 .info1 = sample_info[sample_info$organism == "hyb" &
                                                        sample_info$allele == "cer" &
                                                        sample_info$genotype == "WT",],
                                 .info2 = sample_info[sample_info$organism == "hyb" &
                                                        sample_info$allele == "par" &
                                                        sample_info$genotype == "WT",],
                                 .info3 = sample_info[sample_info$organism == "hyb" &
                                                        sample_info$allele == "cer" &
                                                        sample_info$genotype == paste0(.tf, "delete"),],
                                 .info4 = sample_info[sample_info$organism == "hyb" & 
                                                        sample_info$allele == "par" &
                                                        sample_info$genotype == paste0(.tf, "delete"),],
                                 .normalization = .normalization, .plotlims = .plotlims,
                                 .show_wt = .show_wt, .show_tfdel = .show_tfdel)
  }
  return(p)
}

# Upset plot
# given a group name, which must be a column in .df,
# members of that group (possible values in that column, such as TFs in the deletion column),
# items in the group, other column names in .df (usually gene names, but can be gene-effect-TF combos if organism is group)
# creates an upset plot (3+ group Venn diagram) 
# which enumerates which combinations of groups share how many items
makeUpsetPlot <- function(.df, .group_name, .group_members, .item_names,
                          .min_comb_size = 5) {
  lt <- vector(mode = "list", length = 0)
  for (grpmem in .group_members) {
    lt[[grpmem]] <- filter(.df, .data[[.group_name]] == grpmem)
  }
  lt <- lt |> 
    map(.f = select, .item_names) |> 
    map(.f = \(x) {purrr::reduce(x, .f = paste0)})
  plotdf <- make_comb_mat(lt)
  plotdf <- plotdf[,comb_size(plotdf) >= .min_comb_size]
  p <- UpSet(plotdf, 
             set_order = .group_members,
             comb_order = order(comb_size(plotdf)),
             top_annotation = HeatmapAnnotation( 
               "number of genes" = anno_barplot(comb_size(plotdf), 
                                                ylim = c(0, max(comb_size(plotdf))*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(4, "cm")), 
               annotation_name_side = "left", 
               annotation_name_rot = 90))
  draw(p)
  suppressWarnings(decorate_annotation("number of genes", {
    grid.text(comb_size(plotdf)[column_order(p)], x = seq_along(comb_size(plotdf)), y = unit(comb_size(plotdf)[column_order(p)], "native") + unit(2, "pt"), 
              default.units = "native", just = c("left", "bottom"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
  }))
}

# get a dataframe of gene ontology(goslim) terms found in a given set of genes
getGOSlimDf <- function(.idxs, .group_name, .file_prefix = "gene_ontology/",
                        .min_hits = 5) {
  test_table <- goslim |> filter(ORF %in% .idxs) |> select(GOslim_term) |> table()
  test_table <- test_table[test_table >= .min_hits]
  testdf <- tibble("term" = names(test_table),
                   "group_count" = as.numeric(test_table))
  testdf$overall_count <- map(testdf$term, \(x) {
    ct <- goslim |> filter(GOslim_term == x) |> select(ORF) |> 
      pull() |> unique() |> length()
    return(ct)
  }) |> unlist()
  testdf$exact_pval <- purrr::map2(testdf$group_count, testdf$overall_count, \(x, y) {
    n_mod <- length(.idxs)
    n_genes <- nrow(counts)
    exact_table <- rbind(c(x, n_mod - x), c(y - x, n_genes - n_mod - y))
    exact_result <- fisher.test(exact_table, alternative = "greater")
    return(exact_result$p.value)
  }) |> unlist()
  testdf$sig <- testdf$exact_pval < 0.001
  testdf$genes <- purrr::map(testdf$term, \(x) {
    termgenes <- goslim |> filter(ORF %in% .idxs & GOslim_term == x) |> 
      select(ORF) |> pull() |> unique()
    return(purrr::reduce(termgenes, paste, sep = " "))
  }) |> unlist()
  
  write_csv(arrange(testdf, exact_pval, desc(group_count)), 
            file = paste0(.file_prefix, .group_name, ".csv"),
            quote = "none", 
            col_names = TRUE)
  return(testdf)
}