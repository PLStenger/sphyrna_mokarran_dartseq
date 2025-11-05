#!/bin/bash
#SBATCH --job-name=dartseq_sphyrna
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=500G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/sphyrna_mokarran_dartseq/00_scripts/dartseq_analysis.err"
#SBATCH --output="/home/plstenge/sphyrna_mokarran_dartseq/00_scripts/dartseq_analysis.out"

# Script pour analyses DArTseq sur Sphyrna mokarran
# Date: 2025-11-05
# Auteur: Pierre-Louis Stenger

echo "=========================================="
echo "Démarrage de l'analyse DArTseq"
echo "Date: $(date)"
echo "=========================================="

# Charger le module R (adapter selon votre cluster)
# Vérifier les modules R disponibles avec: module avail R
module load R/4.3.0 || module load R/4.2.0 || module load R/4.1.0 || module load R

# Alternative si R est dans conda
# module load miniconda3
# source activate base

# Vérifier que R est bien chargé
echo "Vérification de R:"
which R
which Rscript
R --version

# Définir les chemins
PROJECT_DIR="/home/plstenge/sphyrna_mokarran_dartseq"
RAW_DATA_DIR="${PROJECT_DIR}/01_raw_data/Report-DSph25-10737"
RESULTS_DIR="${PROJECT_DIR}/03_results"
PLOTS_DIR="${RESULTS_DIR}/plots"
TABLES_DIR="${RESULTS_DIR}/tables"
REPORT_DIR="${RESULTS_DIR}/reports"
SCRIPT_DIR="${PROJECT_DIR}/00_scripts"

# Créer les répertoires de sortie s'ils n'existent pas
mkdir -p ${RESULTS_DIR}
mkdir -p ${PLOTS_DIR}
mkdir -p ${TABLES_DIR}
mkdir -p ${REPORT_DIR}

echo "Répertoires créés:"
echo "  - Résultats: ${RESULTS_DIR}"
echo "  - Plots: ${PLOTS_DIR}"
echo "  - Tables: ${TABLES_DIR}"
echo "  - Rapports: ${REPORT_DIR}"
echo ""

# Lancer le script R
echo "Lancement du script R d'analyse..."
Rscript ${SCRIPT_DIR}/dartseq_analysis.R

# Vérifier si le script R s'est bien exécuté
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Analyse terminée avec succès!"
    echo "Date: $(date)"
    echo "=========================================="
    
    # Créer un rapport HTML si pandoc est disponible
    if command -v pandoc &> /dev/null; then
        echo "Génération du rapport HTML..."
        cd ${REPORT_DIR}
        if [ -f "rapport_dartseq.md" ]; then
            pandoc rapport_dartseq.md -o rapport_dartseq.html --self-contained
            echo "Rapport HTML généré: ${REPORT_DIR}/rapport_dartseq.html"
        fi
    fi
    
    # Créer un résumé des résultats
    SUMMARY_FILE="${REPORT_DIR}/summary.txt"
    echo "Résumé de l'analyse DArTseq - Sphyrna mokarran" > ${SUMMARY_FILE}
    echo "Date: $(date)" >> ${SUMMARY_FILE}
    echo "======================================" >> ${SUMMARY_FILE}
    echo "" >> ${SUMMARY_FILE}
    echo "Fichiers générés:" >> ${SUMMARY_FILE}
    echo "" >> ${SUMMARY_FILE}
    echo "PLOTS:" >> ${SUMMARY_FILE}
    ls -lh ${PLOTS_DIR} >> ${SUMMARY_FILE}
    echo "" >> ${SUMMARY_FILE}
    echo "TABLES:" >> ${SUMMARY_FILE}
    ls -lh ${TABLES_DIR} >> ${SUMMARY_FILE}
    echo "" >> ${SUMMARY_FILE}
    
    # Envoyer le résumé par mail
    if command -v mail &> /dev/null; then
        echo "Envoi du résumé par mail..."
        mail -s "Analyse DArTseq Sphyrna mokarran - Terminée" pierrelouis.stenger@gmail.com < ${SUMMARY_FILE}
    else
        echo "Commande 'mail' non disponible. Résumé sauvegardé dans: ${SUMMARY_FILE}"
    fi
    
    echo ""
    echo "Résultats disponibles dans: ${RESULTS_DIR}"
    
else
    echo ""
    echo "=========================================="
    echo "ERREUR: L'analyse a échoué"
    echo "Date: $(date)"
    echo "Consultez les logs d'erreur: ${SCRIPT_DIR}/dartseq_analysis.err"
    echo "=========================================="
    exit 1
fi
