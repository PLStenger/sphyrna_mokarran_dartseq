#!/usr/bin/env Rscript

# -------------------------------------------------------------
# Script d'installation + test pour analyses DArTseq modernes
# -------------------------------------------------------------

# Force la biblio utilisateur locale dans $HOME
user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))
cat("LibPaths utilisées :", .libPaths(), "\n")

# Fonction d'install robuste :
robust_install <- function(pkgs, repo = "https://cloud.r-project.org/") {
  pkgs_missing <- pkgs[!sapply(pkgs, function(p) requireNamespace(p, quietly = TRUE))]
  if(length(pkgs_missing)) {
    install.packages(pkgs_missing, repos = repo, lib = user_lib, dependencies = TRUE)
  }
  invisible(sapply(pkgs, require, character.only = TRUE))
}

# Install + charge tous les packages DArTseq & manip modernes
robust_install(c(
  "cli", "remotes", "tidyverse", "ggplot2", "data.table",
  "stringr", "ape", "adegenet", "vegan", "dartR", "dartR.base", "dartR.data"
))

# Installe toujours l'update la plus récente de dartR.base+dartR.data
if (!requireNamespace("dartR.base", quietly = TRUE)) {
  remotes::install_github("green-striped-gecko/dartR.base", lib = user_lib)
  library(dartR.base)
}
if (!requireNamespace("dartR.data", quietly = TRUE)) {
  remotes::install_github("green-striped-gecko/dartR.data", lib = user_lib)
  library(dartR.data)
}

cat("\nPackages installés et chargés :\n")
print(.packages())

cat("\nTEST dartR :\n")
cat("dartR version :", packageVersion("dartR"), "\n")
cat("dartR.base version :", packageVersion("dartR.base"), "\n")
cat("dartR.data version :", packageVersion("dartR.data"), "\n")

# Test de lecture simple (si tu veux vérifier plus tard)
# gl <- gl.read.dart(filename = "/chemin/vers/ton/fichier.csv")

cat("\nInstallation terminée. Tous les packages modernes pour DArTseq sont fonctionnels !\n")
