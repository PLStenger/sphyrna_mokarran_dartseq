#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Script R complet pour analyses DArTseq sur Sphyrna mokarran
# Installation robuste des packages et génération automatique des résultats
# Auteur : Pierre-Louis Stenger
# Date   : 2025-11-05
# ------------------------------------------------------------------------------

cat("========================================\n")
cat("Installation des packages nécessaires\n")
cat("========================================\n")

# 1. Installation BiocManager et packages Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}
library(BiocManager)

# Installation des packages Bioconductor requis
bioc_packages <- c("SNPRelate", "gdsfmt")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installation de", pkg, "via Bioconductor...\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# 2. Installation des packages CRAN requis
cran_packages <- c("devtools", "tidyverse", "ggplot2", "stringr", "rmarkdown", 
                   "adegenet", "ade4", "vegan", "ape", "MASS", "reshape2")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installation de", pkg, "via CRAN...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org/", dependencies = TRUE)
  }
}

# 3. Installation spéciale de dartR selon la nouvelle approche dartRverse
cat("Installation de dartRverse...\n")
if (!requireNamespace("dartRverse", quietly = TRUE)) {
  install.packages("dartRverse", repos = "https://cloud.r-project.org/")
}

# Chargement des librairies essentielles
suppressMessages({
  library(dartRverse)
  library(tidyverse)
  library(ggplot2)
  library(stringr)
})

cat("Packages installés avec succès!\n\n")

# 4. Définition des chemins
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

cat("========================================\n")
cat("Import et analyse des données DArTseq\n")
cat("========================================\n")

# 5. Import des données DArTseq
setwd(raw_data_dir)

# Vérifier la présence du fichier
snp_file <- "Report_DSph25-10737_2_moreOrders_SNP_2.csv"
if (!file.exists(snp_file)) {
  stop("Fichier SNP non trouvé : ", snp_file)
}

cat("Import du fichier :", snp_file, "\n")
gl <- gl.read.dart(filename = snp_file, ind.metafile = NULL)

# Informations de base
n_ind <- nInd(gl)
n_loc <- nLoc(gl)
n_pop <- nPop(gl)

cat("Nombre d'individus    :", n_ind, "\n")
cat("Nombre de loci        :", n_loc, "\n")
cat("Nombre de populations :", n_pop, "\n\n")

# 6. Sauvegarde des informations de base
basic_info <- data.frame(
  Metric = c("Nombre_individus", "Nombre_loci", "Nombre_populations"),
  Value = c(n_ind, n_loc, n_pop)
)
write.csv(basic_info, file.path(tables_dir, "basic_info.csv"), row.names = FALSE)

# 7. Nettoyage et correction des noms (adapté de votre script original)
ind_names <- indNames(gl)

# Dictionnaire de noms (exemple - à adapter selon vos données)
name_dict <- c(
  "MOKARR" = "MOKARR",
  "HAMMERHEAD" = "HAMMER",
  "SPHYR" = "SPHYR"
  # Ajoutez vos corrections ici
)

# Application des corrections si nécessaire
corrected_names <- ind_names
for (old_name in names(name_dict)) {
  corrected_names <- str_replace_all(corrected_names, old_name, name_dict[old_name])
}
indNames(gl) <- corrected_names

# Définition des populations basée sur les préfixes (exemple)
if (n_pop <= 1) {
  # Créer des populations basées sur les premiers caractères des noms
  pop_codes <- str_sub(corrected_names, 1, 3)
  pop(gl) <- as.factor(pop_codes)
  cat("Populations redéfinies basées sur les préfixes des noms\n")
}

# 8. Filtrages qualité
cat("========================================\n")
cat("Filtrages qualité\n")
cat("========================================\n")

# Statistiques avant filtrage
cat("Avant filtrage - Loci:", nLoc(gl), "Individus:", nInd(gl), "\n")

# Filtrage par RepAvg (si disponible)
if ("RepAvg" %in% names(gl@other$loc.metrics)) {
  gl <- gl.filter.repavg(gl, threshold = 0.95, verbose = TRUE)
  cat("Après filtrage RepAvg - Loci:", nLoc(gl), "\n")
}

# Filtrage par CallRate
gl <- gl.filter.callrate(gl, method = "loc", threshold = 0.80, verbose = TRUE)
cat("Après filtrage CallRate - Loci:", nLoc(gl), "\n")

# Suppression des loci monomorphes
gl <- gl.filter.monomorphs(gl, verbose = TRUE)
cat("Après suppression monomorphes - Loci:", nLoc(gl), "\n")

# Sauvegarde des statistiques de filtrage
filtering_stats <- data.frame(
  Step = c("Initial", "After_RepAvg", "After_CallRate", "After_Monomorphs"),
  N_Loci = c(n_loc, nLoc(gl), nLoc(gl), nLoc(gl)),
  N_Individuals = c(n_ind, nInd(gl), nInd(gl), nInd(gl))
)
write.csv(filtering_stats, file.path(tables_dir, "filtering_statistics.csv"), row.names = FALSE)

# 9. Analyses de génétique des populations
cat("========================================\n")
cat("Analyses génétiques\n")
cat("========================================\n")

# 9.1 Calcul de l'hétérozygotie
cat("Calcul de l'hétérozygotie...\n")
het_results <- gl.report.heterozygosity(gl, method = "pop")
write.csv(het_results, file.path(tables_dir, "heterozygosity_by_pop.csv"), row.names = FALSE)

# 9.2 PCoA
cat("Analyse en coordonnées principales (PCoA)...\n")
pcoa_result <- gl.pcoa(gl, nf = 3)

# Extraction et sauvegarde des scores PCoA
pcoa_scores <- as.data.frame(pcoa_result$scores)
pcoa_scores$Individual <- rownames(pcoa_scores)
pcoa_scores$Population <- pop(gl)
write.csv(pcoa_scores, file.path(tables_dir, "pcoa_scores.csv"), row.names = FALSE)

# 9.3 Calcul de Fst entre populations (si plus d'une population)
if (nPop(gl) > 1) {
  cat("Calcul des Fst entre populations...\n")
  fst_matrix <- gl.fst.pop(gl)
  write.csv(as.matrix(fst_matrix), file.path(tables_dir, "fst_matrix.csv"), row.names = TRUE)
}

# 10. Génération des graphiques
cat("========================================\n")
cat("Génération des graphiques\n")
cat("========================================\n")

# 10.1 Plot PCoA
p_pcoa <- ggplot(pcoa_scores, aes(x = Axis1, y = Axis2, color = Population)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(aes(group = Population), type = "norm", level = 0.95) +
  labs(title = "Analyse en Coordonnées Principales (PCoA)",
       x = paste0("PC1 (", round(pcoa_result$eig[1]/sum(pcoa_result$eig)*100, 1), "%)"),
       y = paste0("PC2 (", round(pcoa_result$eig[2]/sum(pcoa_result$eig)*100, 1), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(plots_dir, "pcoa_analysis.png"), p_pcoa, width = 10, height = 8, dpi = 300)

# 10.2 Barplot hétérozygotie
if (nrow(het_results) > 1) {
  p_het <- ggplot(het_results, aes(x = pop, y = Ho)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(ymin = Ho - 0.01, ymax = Ho + 0.01), width = 0.2) +
    labs(title = "Hétérozygotie Observée par Population",
         x = "Population", y = "Hétérozygotie Observée (Ho)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(plots_dir, "heterozygosity_barplot.png"), p_het, width = 8, height = 6, dpi = 300)
}

# 10.3 Heatmap Fst (si plus d'une population)
if (nPop(gl) > 1 && exists("fst_matrix")) {
  fst_df <- expand.grid(Pop1 = rownames(fst_matrix), Pop2 = colnames(fst_matrix))
  fst_df$Fst <- as.vector(as.matrix(fst_matrix))
  
  p_fst <- ggplot(fst_df, aes(x = Pop1, y = Pop2, fill = Fst)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = "Matrice de Fst entre Populations", x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(plots_dir, "fst_heatmap.png"), p_fst, width = 8, height = 6, dpi = 300)
}

# 11. Création du rapport
cat("========================================\n")
cat("Génération du rapport\n")
cat("========================================\n")

# Création du rapport markdown
report_file <- file.path(report_dir, "rapport_dartseq.md")

# En-tête du rapport
cat("# Rapport d'Analyse DArTseq - Sphyrna mokarran\n\n", file = report_file)
cat("**Date d'analyse :** ", as.character(Sys.Date()), "\n", file = report_file, append = TRUE)
cat("**Auteur :** Pierre-Louis Stenger\n\n", file = report_file, append = TRUE)

# Résumé des données
cat("## Résumé des données\n\n", file = report_file, append = TRUE)
cat("- **Nombre d'individus :** ", nInd(gl), "\n", file = report_file, append = TRUE)
cat("- **Nombre de loci (après filtrage) :** ", nLoc(gl), "\n", file = report_file, append = TRUE)
cat("- **Nombre de populations :** ", nPop(gl), "\n\n", file = report_file, append = TRUE)

# Liste des fichiers générés
cat("## Fichiers générés\n\n", file = report_file, append = TRUE)
cat("### Tables (.csv)\n", file = report_file, append = TRUE)
table_files <- list.files(tables_dir, pattern = "\\.csv$")
for (f in table_files) {
  cat("- ", f, "\n", file = report_file, append = TRUE)
}

cat("\n### Graphiques (.png)\n", file = report_file, append = TRUE)
plot_files <- list.files(plots_dir, pattern = "\\.png$")
for (f in plot_files) {
  cat("- ", f, "\n", file = report_file, append = TRUE)
}

# Conversion en HTML si possible
html_file <- file.path(report_dir, "rapport_dartseq.html")
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  try({
    rmarkdown::render(report_file, output_file = html_file, quiet = TRUE)
    cat("Rapport HTML généré :", html_file, "\n")
  }, silent = TRUE)
}

# Résumé final
cat("\n========================================\n")
cat("ANALYSE TERMINÉE AVEC SUCCÈS!\n")
cat("========================================\n")
cat("Résultats disponibles dans :", results_dir, "\n")
cat("- Tables :", length(list.files(tables_dir)), "fichiers\n")
cat("- Graphiques :", length(list.files(plots_dir)), "fichiers\n")
cat("- Rapport :", report_file, "\n")

# Sauvegarde de l'objet genlight final
save(gl, file = file.path(results_dir, "genlight_filtered.RData"))
cat("Objet genlight sauvegardé : genlight_filtered.RData\n")
