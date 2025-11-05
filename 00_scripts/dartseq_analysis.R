#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Script R complet pour analyses DArTseq sur Sphyrna mokarran
# Installation dans répertoire utilisateur avec gestion des permissions
# Auteur : Pierre-Louis Stenger
# Date   : 2025-11-05
# ------------------------------------------------------------------------------

cat("========================================\n")
cat("Configuration de l'environnement R\n")
cat("========================================\n")

# 1. Définir une bibliothèque utilisateur pour les packages
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
}

# Créer le répertoire s'il n'existe pas
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
  cat("Bibliothèque utilisateur créée:", user_lib, "\n")
}

# Ajouter au chemin de recherche des bibliothèques
.libPaths(c(user_lib, .libPaths()))
cat("Bibliothèques R utilisées:\n")
print(.libPaths())

cat("\n========================================\n")
cat("Installation des packages nécessaires\n")
cat("========================================\n")

# 2. Fonction d'installation sécurisée
install_if_missing <- function(pkg, repo = "CRAN") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installation de", pkg, "depuis", repo, "...\n")
    if (repo == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org/")
      }
      BiocManager::install(pkg, lib = user_lib, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org/", dependencies = TRUE)
    }
  } else {
    cat(pkg, "déjà installé\n")
  }
}

# 3. Installation des packages Bioconductor
cat("\n--- Installation packages Bioconductor ---\n")
install_if_missing("BiocManager", "CRAN")
library(BiocManager)
install_if_missing("SNPRelate", "Bioconductor")
install_if_missing("gdsfmt", "Bioconductor")

# 4. Installation des packages CRAN de base
cat("\n--- Installation packages CRAN ---\n")
essential_packages <- c("tidyverse", "ggplot2", "stringr", "data.table", 
                        "ade4", "vegan", "ape", "MASS", "reshape2")
for (pkg in essential_packages) {
  install_if_missing(pkg, "CRAN")
}

# 5. Installation de adegenet (dépendance importante)
cat("\n--- Installation adegenet ---\n")
install_if_missing("adegenet", "CRAN")

# 6. Installation de dartR.data et dartR.base
cat("\n--- Installation dartR ---\n")
install_if_missing("dartR.data", "CRAN")
install_if_missing("dartR.base", "CRAN")

# 7. Chargement des librairies
cat("\n========================================\n")
cat("Chargement des librairies\n")
cat("========================================\n")

suppressPackageStartupMessages({
  library(dartR.base)
  library(tidyverse)
  library(ggplot2)
  library(stringr)
})

cat("Packages chargés avec succès!\n")

# 8. Définition des chemins
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

cat("\n========================================\n")
cat("Import et analyse des données DArTseq\n")
cat("========================================\n")

# 9. Import des données DArTseq
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

# 10. Sauvegarde des informations de base
basic_info <- data.frame(
  Metric = c("Nombre_individus", "Nombre_loci", "Nombre_populations"),
  Value = c(n_ind, n_loc, n_pop)
)
write.csv(basic_info, file.path(tables_dir, "basic_info.csv"), row.names = FALSE)

# 11. Nettoyage et correction des noms d'individus
cat("Correction des noms d'individus...\n")
ind_names <- indNames(gl)

# Dictionnaire de correction (à adapter selon vos besoins)
# Exemple simple : nettoyer les espaces et caractères spéciaux
ind_names_clean <- str_replace_all(ind_names, "[^[:alnum:]_-]", "_")
indNames(gl) <- ind_names_clean

# Définition des populations si besoin
if (n_pop <= 1) {
  # Créer des populations basées sur les 3 premiers caractères
  pop_codes <- str_sub(ind_names_clean, 1, 3)
  pop(gl) <- as.factor(pop_codes)
  cat("Populations définies basées sur les préfixes des noms\n")
  cat("Populations identifiées:", length(unique(pop(gl))), "\n")
}

# 12. Filtrages qualité
cat("\n========================================\n")
cat("Filtrages qualité\n")
cat("========================================\n")

cat("Avant filtrage - Loci:", nLoc(gl), "Individus:", nInd(gl), "\n")

# Filtrage par RepAvg si disponible
if ("RepAvg" %in% names(gl@other$loc.metrics)) {
  cat("Filtrage RepAvg >= 0.95...\n")
  gl <- gl.filter.repavg(gl, threshold = 0.95, verbose = 3)
  cat("Après filtrage RepAvg - Loci:", nLoc(gl), "\n")
}

# Filtrage par CallRate
cat("Filtrage CallRate >= 0.80...\n")
gl <- gl.filter.callrate(gl, method = "loc", threshold = 0.80, verbose = 3)
cat("Après filtrage CallRate - Loci:", nLoc(gl), "\n")

# Suppression des loci monomorphes
cat("Suppression des loci monomorphes...\n")
gl <- gl.filter.monomorphs(gl, verbose = 3)
cat("Après suppression monomorphes - Loci:", nLoc(gl), "\n")

# Sauvegarde des statistiques de filtrage
filtering_stats <- data.frame(
  Step = c("Initial", "Final"),
  N_Loci = c(n_loc, nLoc(gl)),
  N_Individuals = c(n_ind, nInd(gl))
)
write.csv(filtering_stats, file.path(tables_dir, "filtering_statistics.csv"), row.names = FALSE)

# 13. Analyses génétiques
cat("\n========================================\n")
cat("Analyses de génétique des populations\n")
cat("========================================\n")

# 13.1 Hétérozygotie
cat("Calcul de l'hétérozygotie...\n")
het_results <- gl.report.heterozygosity(gl, method = "pop")
write.csv(het_results, file.path(tables_dir, "heterozygosity_by_pop.csv"), row.names = FALSE)

# 13.2 PCoA
cat("Analyse en Coordonnées Principales (PCoA)...\n")
pcoa_result <- gl.pcoa(gl, nfactors = 3)

# Extraction des scores
pcoa_scores <- as.data.frame(pcoa_result$scores)
pcoa_scores$Individual <- rownames(pcoa_scores)
pcoa_scores$Population <- as.character(pop(gl))
write.csv(pcoa_scores, file.path(tables_dir, "pcoa_scores.csv"), row.names = FALSE)

# 13.3 Fst entre populations (si applicable)
if (nPop(gl) > 1) {
  cat("Calcul des Fst entre populations...\n")
  fst_matrix <- gl.fst.pop(gl)
  write.csv(as.matrix(fst_matrix), file.path(tables_dir, "fst_matrix.csv"), row.names = TRUE)
}

# 14. Génération des graphiques
cat("\n========================================\n")
cat("Génération des graphiques\n")
cat("========================================\n")

# 14.1 PCoA plot
var_explained <- pcoa_result$eig / sum(pcoa_result$eig) * 100

p_pcoa <- ggplot(pcoa_scores, aes(x = Axis1, y = Axis2, color = Population)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "Analyse en Coordonnées Principales (PCoA)",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(plots_dir, "pcoa_analysis.png"), p_pcoa, width = 10, height = 8, dpi = 300)
cat("PCoA plot sauvegardé\n")

# 14.2 Hétérozygotie barplot
if (nrow(het_results) > 1) {
  p_het <- ggplot(het_results, aes(x = pop, y = Ho)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    labs(title = "Hétérozygotie Observée par Population",
         x = "Population", y = "Hétérozygotie Observée (Ho)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(plots_dir, "heterozygosity_barplot.png"), p_het, width = 8, height = 6, dpi = 300)
  cat("Hétérozygotie barplot sauvegardé\n")
}

# 14.3 Heatmap Fst
if (nPop(gl) > 1 && exists("fst_matrix")) {
  fst_df <- expand.grid(Pop1 = rownames(fst_matrix), Pop2 = colnames(fst_matrix))
  fst_df$Fst <- as.vector(as.matrix(fst_matrix))
  
  p_fst <- ggplot(fst_df, aes(x = Pop1, y = Pop2, fill = Fst)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = median(fst_df$Fst, na.rm = TRUE)) +
    labs(title = "Matrice de Fst entre Populations", x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(plots_dir, "fst_heatmap.png"), p_fst, width = 8, height = 6, dpi = 300)
  cat("Fst heatmap sauvegardé\n")
}

# 15. Création du rapport
cat("\n========================================\n")
cat("Génération du rapport\n")
cat("========================================\n")

report_file <- file.path(report_dir, "rapport_dartseq.md")

cat("# Rapport d'Analyse DArTseq - Sphyrna mokarran\n\n", file = report_file)
cat("**Date d'analyse :** ", as.character(Sys.Date()), "\n", file = report_file, append = TRUE)
cat("**Auteur :** Pierre-Louis Stenger\n\n", file = report_file, append = TRUE)

cat("## Résumé des données\n\n", file = report_file, append = TRUE)
cat("- **Nombre d'individus :** ", nInd(gl), "\n", file = report_file, append = TRUE)
cat("- **Nombre de loci (après filtrage) :** ", nLoc(gl), "\n", file = report_file, append = TRUE)
cat("- **Nombre de populations :** ", nPop(gl), "\n\n", file = report_file, append = TRUE)

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

# Sauvegarde de l'objet genlight final
save(gl, file = file.path(results_dir, "genlight_filtered.RData"))

cat("\n========================================\n")
cat("ANALYSE TERMINÉE AVEC SUCCÈS!\n")
cat("========================================\n")
cat("Résultats disponibles dans :", results_dir, "\n")
cat("- Tables :", length(list.files(tables_dir)), "fichiers\n")
cat("- Graphiques :", length(list.files(plots_dir)), "fichiers\n")
cat("- Rapport :", report_file, "\n")
cat("- Objet genlight : genlight_filtered.RData\n")
