# Load packages
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(stringi)

# Set seed
set.seed("1234")

# GSEA function
GSEA <- function(gene_file, geneSet_name, GO_file, pval = 0.01, collapsing = TRUE) {
  genes <- read_table(gene_file)
  genes <- arrange(genes, ZSTAT)
  gene_list <- genes$ZSTAT
  names(gene_list) <- genes$SYMBOL
  myGO <- fgsea::gmtPathways(GO_file)

  fgRes <- fgsea::fgsea(
    pathways = myGO,
    stats = gene_list,
    minSize = 15, ## minimum gene set size
    maxSize = 500
  ) %>%
    as.data.frame() %>%
    dplyr::filter(pval < !!pval) %>%
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  if (collapsing) {
    message("Collapsing Pathways -----")
    concise_pathways <- collapsePathways(data.table::as.data.table(fgRes),
      pathways = myGO,
      stats = gene_list
    )
    fgRes_collapsed <- fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
    if (nrow(fgRes_collapsed) == 0) {
      message("No collapsable terms found.")
    } else {
      fgRes <- fgRes_collapsed
      message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
    }
  }
  fgRes$Enrichment <- ifelse(fgRes$NES > 0, "Pos. Z-score", "Neg. Z-score")
  fgRes$GeneSet <- geneSet_name
  filtRes <- rbind(
    head(fgRes, n = 10),
    tail(fgRes, n = 10)
  )

  total_up <- sum(fgRes$Enrichment == "Pos. Z-score")
  total_down <- sum(fgRes$Enrichment == "Neg. Z-score")
  header <- paste0(geneSet_name, " - top pathways (total: up=", total_up, ", down=", total_down, ")")

  colos <- setNames(
    c("firebrick", "#112f43"),
    c("Pos. Z-score", "Neg. Z-score")
  )

  filtRes$pathway <- stri_replace_all_regex(filtRes$pathway, "_", " ")

  g1 <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(colour = Enrichment, size = size)) +
    scale_colour_manual(values = colos) +
    geom_hline(yintercept = 0) +
    scale_size_continuous(range = c(0.5, 3)) +
    coord_flip() +
    labs(
      x = "Pathway", y = "Normalized enrichment score",
      title = header
    ) +
    theme_bw(base_size = 7, base_family = "sans") +
    theme(
      plot.title = element_text(size = 7),
      axis.text = element_text(size = 5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 5),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5),
      legend.key.height = unit(10, "pt"),
      panel.grid = element_line(linewidth = .25)
    ) +
    guides(size = guide_legend(title = "Size"))

  output <- list("Results" = fgRes, "Plot" = g1)
  return(output)
}
