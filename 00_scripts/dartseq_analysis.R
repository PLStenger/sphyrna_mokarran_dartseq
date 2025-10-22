#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Script R complet pour analyses DArTseq sur Sphyrna mokarran
# Génération automatique des tables, plots et rapport Markdown
# Auteur : Pierre-Louis Stenger
# Date   : 2025-10-22
# ------------------------------------------------------------------------------

# 1. Installation et chargement des packages
packages <- c("dartR", "tidyverse", "ggplot2", "stringr", "rmarkdown")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# 2. Définition des chemins
project_dir   <- "/home/plstenge/sphyrna_mokarran_dartseq"
raw_data_dir  <- file.path(project_dir, "01_raw_data/Report-DSph25-10737")
results_dir   <- file.path(project_dir, "03_results")
plots_dir     <- file.path(results_dir, "plots")
tables_dir    <- file.path(results_dir, "tables")
report_dir    <- file.path(results_dir, "reports")

# Création des répertoires de sortie
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir,  showWarnings = FALSE, recursive = TRUE)

# 3. Import des données DArTseq
setwd(raw_data_dir)
gl <- gl.read.dart(
  filename     = "Report_DSph25-10737_2_moreOrders_SNP_2.csv",
  ind.metafile = NULL
)
cat("\nNombre d'individus :", nInd(gl), "\n")
cat("Nombre de loci        :", nLoc(gl), "\n")
cat("Nombre de populations :", nPop(gl), "\n\n")

# 4. Correction des noms d’individus et définition des populations
# (Adapter ce bloc selon ton dictionnaire de noms)
ind_names <- indNames(gl)
# Exemple : ind_names <- str_replace(ind_names, "OLD", "NEW")
# indNames(gl) <- ind_names
# pop(gl)      <- factor(str_sub(ind_names, 1, 3))

# 5. Filtrages qualité
# Filtrer par moyenne de répétabilité (RepAvg) >= 0.95
gl <- gl.filter.repavg(gl, threshold = 0.95, verbose = TRUE)

# Filtrer par taux d'appel (CallRate) >= 0.80
gl <- gl.filter.callrate(gl, threshold = 0.80, verbose = TRUE)

# Enlever les loci monomorphes
gl <- gl.filter.monomorphs(gl)

# 6. Analyses population
# 6.1 PCoA
pcoa_res <- gl.pcoa(gl, method = "euclidean", nf = 3)
# Export table PCoA
write.csv(pcoa_res$scores,
          file.path(tables_dir, "pcoa_scores.csv"),
          row.names = TRUE)

# 6.2 Fst par paire de populations
fst_mat <- gl.fst.pop(gl)
write.csv(fst_mat,
          file.path(tables_dir, "fst_matrix.csv"),
          row.names = TRUE)

# 6.3 Hétérozygotie
het <- gl.report.heterozygosity(gl)
write.csv(het,
          file.path(tables_dir, "heterozygosity.csv"),
          row.names = FALSE)

# 6.4 Test de Hardy–Weinberg par locus
hwe <- gl.filter.hwe(gl, alpha = 0.05, verbose = TRUE)
write.csv(hwe,
          file.path(tables_dir, "hwe_test.csv"),
          row.names = FALSE)

# 7. Génération des plots
# 7.1 PCoA plot
pcoa_df <- as_tibble(pcoa_res$scores) %>%
  mutate(Ind = rownames(pcoa_res$scores))
p <- ggplot(pcoa_df, aes(PCoA1, PCoA2, label = Ind)) +
  geom_point(size = 3) +
  geom_text_repel() +
  theme_minimal() +
  labs(title = "PCoA des individus", x = "PCoA1", y = "PCoA2")
ggsave(file.path(plots_dir, "pcoa_plot.png"), plot = p, width = 8, height = 6)

# 7.2 Heatmap Fst
fst_melt <- as.data.frame(as.table(fst_mat))
colnames(fst_melt) <- c("Pop1", "Pop2", "Fst")
h <- ggplot(fst_melt, aes(Pop1, Pop2, fill = Fst)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Matrix Fst", x = "", y = "")
ggsave(file.path(plots_dir, "fst_heatmap.png"), plot = h, width = 6, height = 5)

# 8. Fonction utilitaire pour dumping de tables
dump_table <- function(obj, name) {
  write.csv(obj,
            file.path(tables_dir, paste0(name, ".csv")),
            row.names = FALSE)
}
# Exemples
# dump_table(het,        "heterozygosity_results")
# dump_table(fst_mat,    "fst_matrix_full")

# 9. Création du rapport Markdown
report_md <- file.path(report_dir, "rapport_dartseq.md")
cat("# Rapport Analyse DArTseq\n\n", file = report_md)
cat("**Date** : ", Sys.Date(), "\n\n", file = report_md, append = TRUE)
cat("## Résumé\n", file = report_md, append = TRUE)
cat("Analyse SNPs sur Sphyrna mokarran réalisée via dartR.\n\n",
    file = report_md, append = TRUE)
cat("## Fichiers générés\n\n", file = report_md, append = TRUE)
cat("### Tables\n", file = report_md, append = TRUE)
cat(paste0("- ", list.files(tables_dir), collapse = "\n"),
    "\n\n", file = report_md, append = TRUE)
cat("### Plots\n", file = report_md, append = TRUE)
cat(paste0("- ", list.files(plots_dir), collapse = "\n"),
    "\n", file = report_md, append = TRUE)

# 10. Optionnel : conversion en HTML et envoi par mail
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  try({
    rmarkdown::render(input  = report_md,
                      output_file = file.path(report_dir, "rapport_dartseq.html"),
                      quiet = TRUE)
  }, silent = TRUE)
}

# Envoi du rapport par mail si la commande mail est dispo
if (Sys.info()["sysname"] == "Linux" && nzchar(Sys.which("mail"))) {
  system2("mail",
          args = c("-s", shQuote("Rapport DArTseq Sphyrna mokarran"),
                   "pierrelouis.stenger@gmail.com"),
          stdin = report_md)
}

cat("\n=== Analyse terminée ===\n")
cat("Résultats dans :", results_dir, "\n")
