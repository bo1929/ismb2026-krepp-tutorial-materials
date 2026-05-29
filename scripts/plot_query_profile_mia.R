#!/usr/bin/env Rscript
#
# Visualize Kraken-style profile.tsv for ONE sample across taxonomy ranks.
# Uses mia::agglomerateByRank + stacked bars (same compositional view as
# miaViz::plotAbundance). Multiple ranks appear as ggplot facets for a compact
# layout (see https://microbiome.github.io/miaViz/reference/plotAbundance.html).
#
# Usage (from repo root):
#   Rscript scripts/plot_query_profile_mia.R [profile.tsv] [out.png]
#

args <- commandArgs(trailingOnly = TRUE)
profile_path <- if (length(args) >= 1) args[[1]] else "data-new/profile.tsv"
out_path <- if (length(args) >= 2) args[[2]] else "figures/profile_mixture_taxonomy.png"

.libPaths(c(normalizePath(".R/library", mustWork = FALSE), .libPaths()))

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(mia)
  library(ggplot2)
  library(colorspace)
  library(dplyr)
})

rank_slots <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

parse_profile_species <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hdr_idx <- grep("^@@TAXID", lines)
  if (!length(hdr_idx)) {
    stop("Missing @@TAXID header in profile.tsv")
  }
  hdr <- sub("^@@", "", lines[[hdr_idx[[1]]]])
  cn <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
  body <- lines[-seq_len(hdr_idx[[1]])]
  body <- body[nchar(trimws(body)) > 0]
  rows <- strsplit(body, "\t", fixed = TRUE)
  lens <- vapply(rows, length, 1L)
  rows <- rows[lens >= length(cn)]
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- cn[seq_len(ncol(df))]
  df <- df[df$RANK == "species", ]
  df$PERCENTAGE <- as.numeric(df$PERCENTAGE)
  df
}

strip_numeric_tokens <- function(tok) {
  tok <- trimws(tok)
  tok <- tok[nchar(tok) > 0]
  tok[!grepl("^[0-9]+$", tok)]
}

tokens_from_taxpathsn <- function(taxpathsn) {
  strip_numeric_tokens(strsplit(taxpathsn, "|", fixed = TRUE)[[1]])
}

assign_ranks <- function(tok) {
  out <- stats::setNames(rep(NA_character_, length(rank_slots)), rank_slots)
  n <- length(tok)
  if (n == 0) {
    return(as.list(out))
  }
  if (n >= 7) {
    for (i in seq_len(7)) {
      out[[rank_slots[i]]] <- tok[[i]]
    }
  } else if (n >= 2) {
    out[["Kingdom"]] <- tok[[1]]
    out[["Species"]] <- tok[[n]]
    out[["Genus"]] <- tok[[n - 1]]
    if (n >= 4) {
      out[["Phylum"]] <- tok[[2]]
      out[["Class"]] <- tok[[3]]
    }
    if (n >= 6) {
      out[["Order"]] <- tok[[4]]
      out[["Family"]] <- tok[[5]]
    }
  }
  as.list(out)
}

species_to_tse <- function(df) {
  tax_rows <- lapply(df$TAXPATHSN, function(s) assign_ranks(tokens_from_taxpathsn(s)))
  rd <- do.call(rbind, lapply(tax_rows, function(x) data.frame(x, stringsAsFactors = FALSE)))
  rownames(rd) <- paste0("sp_", df$TAXID)
  assay_mat <- matrix(df$PERCENTAGE, ncol = 1, dimnames = list(rownames(rd), "mixture"))
  SummarizedExperiment(assays = list(counts = assay_mat), rowData = DataFrame(rd))
}

build_phylum_palette <- function(rd_df) {
  phy <- sort(unique(rd_df$Phylum[!is.na(rd_df$Phylum) & nchar(rd_df$Phylum) > 0]))
  if (!length(phy)) {
    return(character(0))
  }
  cols <- qualitative_hcl(length(phy), palette = "Dark 3")
  stats::setNames(cols, phy)
}

species_colors <- function(rd_df, abund_vec, phylum_palette) {
  n <- nrow(rd_df)
  out <- rep("#b4b4b4", n)
  names(out) <- rownames(rd_df)
  ph_groups <- split(seq_len(n), rd_df$Phylum)
  for (phy in names(ph_groups)) {
    idx <- ph_groups[[phy]]
    base <- phylum_palette[[phy]]
    if (is.null(base)) {
      base <- "#666666"
    }
    m <- length(idx)
    idx <- idx[order(abund_vec[idx], decreasing = TRUE)]
    if (m == 1) {
      out[idx] <- base
    } else {
      cr <- coords(as(hex2RGB(base), "polarLUV"))
      out[idx] <- sequential_hcl(
        m,
        h1 = cr["H"],
        c1 = max(26, cr["C"] * 0.92),
        l1 = cr["L"],
        l2 = min(93, cr["L"] + 32)
      )
    }
  }
  out
}

weighted_hex_blend <- function(hex_vec, w_vec) {
  if (!length(hex_vec)) {
    return("#c8c8c8")
  }
  w_vec <- as.numeric(w_vec)
  ok <- is.finite(w_vec) & w_vec > 0 & !is.na(hex_vec)
  hex_vec <- hex_vec[ok]
  w_vec <- w_vec[ok]
  if (!length(hex_vec)) {
    return("#c8c8c8")
  }
  w_vec <- w_vec / sum(w_vec)
  rgb_mat <- t(col2rgb(hex_vec)) * w_vec
  vals <- colSums(rgb_mat)
  grDevices::rgb(vals[1], vals[2], vals[3], maxColorValue = 255)
}

taxon_color_at_rank <- function(rd_df, abund_vec, species_hex, rk, taxon_label) {
  if (is.na(taxon_label) || !nzchar(taxon_label)) {
    return("#d9d9d9")
  }
  idx <- which(rd_df[[rk]] == taxon_label & !is.na(rd_df[[rk]]))
  if (!length(idx)) {
    return("#d9d9d9")
  }
  weighted_hex_blend(species_hex[idx], abund_vec[idx])
}

rank_plot_frame <- function(tse, rk, rd_df, abund_vec, species_hex) {
  tx <- agglomerateByRank(tse, rank = rk)
  if (!nrow(tx)) {
    return(NULL)
  }
  tx <- transformAssay(
    tx,
    assay.type = "counts",
    method = "relabundance",
    name = "rel"
  )
  y <- as.numeric(assay(tx, "rel")[, 1])
  lab <- as.character(rowData(tx)[[rk]])
  ok <- !is.na(lab) & nzchar(lab) & y > 0
  data.frame(
    sample = colnames(tx)[[1]],
    rank = rk,
    taxon = lab[ok],
    abundance = y[ok],
    stringsAsFactors = FALSE
  )
}

main <- function() {
  df <- parse_profile_species(profile_path)
  if (!nrow(df)) {
    stop("No species rows in ", profile_path)
  }

  tse <- species_to_tse(df)
  rownames(tse) <- make.unique(rownames(tse))

  rd_df <- as.data.frame(rowData(tse))
  rownames(rd_df) <- rownames(tse)
  abund_vec <- as.numeric(assay(tse, "counts")[, 1])
  names(abund_vec) <- rownames(tse)

  phylum_palette <- build_phylum_palette(rd_df)
  sp_hex <- species_colors(rd_df, abund_vec, phylum_palette)
  kingdom_cols <- c(Archaea = "#9e4059", Bacteria = "#366391")

  ranks_plot <- intersect(taxonomyRanks(tse), rank_slots)

  frames <- lapply(ranks_plot, function(rk) {
    rank_plot_frame(tse, rk, rd_df, abund_vec, sp_hex)
  })
  plot_df <- bind_rows(frames)
  if (!nrow(plot_df)) {
    stop("No abundance rows to plot.")
  }

  plot_df$rank <- factor(plot_df$rank, levels = ranks_plot)

  plot_df <- plot_df |>
    group_by(rank) |>
    mutate(abundance = abundance / sum(abundance)) |>
    ungroup()

  plot_df <- plot_df |>
    mutate(
      fill_key = paste(as.character(rank), taxon, sep = " | "),
      fill_col = vapply(seq_len(nrow(plot_df)), function(i) {
        rk <- as.character(plot_df$rank[[i]])
        tx <- plot_df$taxon[[i]]
        if (rk == "Kingdom" && tx %in% names(kingdom_cols)) {
          return(unname(kingdom_cols[[tx]]))
        }
        taxon_color_at_rank(rd_df, abund_vec, sp_hex, rk, tx)
      }, character(1))
    )

  fill_tbl <- plot_df |> distinct(fill_key, fill_col)
  fill_map <- stats::setNames(fill_tbl$fill_col, fill_tbl$fill_key)

  plot_df <- plot_df |>
    group_by(rank) |>
    arrange(desc(abundance), taxon) |>
    mutate(taxon = factor(taxon, levels = unique(taxon))) |>
    ungroup() |>
    mutate(fill_key = factor(fill_key, levels = names(fill_map)))

  n_facets <- nlevels(plot_df$rank)
  ncol_facet <- min(4L, n_facets)

  p <- ggplot(plot_df, aes(x = sample, y = abundance, fill = fill_key)) +
    geom_col(width = 0.85, color = NA) +
    facet_wrap(vars(rank), ncol = ncol_facet, scales = "fixed") +
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(
      values = fill_map,
      breaks = names(fill_map),
      labels = sub("^[^|]+\\| ", "", names(fill_map)),
      name = NULL
    ) +
    labs(
      title = "Simulated mixture (one sample): relative abundance by rank",
      x = NULL,
      y = "Relative abundance"
    ) +
    theme_bw(base_size = 9) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      strip.background = element_rect(fill = "grey92", color = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 7.5),
      legend.key.height = unit(0.32, "cm"),
      legend.key.width = unit(0.4, "cm"),
      plot.title = element_text(face = "bold", size = 11),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 0.65)))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  h <- max(3.1, 1.22 * ceiling(n_facets / ncol_facet))
  ggsave(out_path, p, width = 8.6, height = h, dpi = 160, limitsize = FALSE)
  message("Wrote ", normalizePath(out_path, mustWork = FALSE))
}

main()
