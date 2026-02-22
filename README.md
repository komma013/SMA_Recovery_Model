# SMA Sequential Therapy Recovery Model

This repository contains the R code used to generate the model predictions and figures presented in the science fair poster:

These scripts implement a biologically constrained recovery model to estimate cohort-level trajectories of motor scores following sequential SMA therapies.

Predicted trajectories are generated for:

* **CHOP-INTEND** — motor function scale for Type I SMA
* **HFMSE** — motor function scale for Type II SMA

Predictions are evaluated at 6, 12, 18, and 24 months after initiation of the second therapy.

*Files:

* `CHOP_Model.R` — implements the CHOP-INTEND recovery model and generates the CHOP trajectory figures shown on the poster
* `HFMSE_Model.R` — implements the HFMSE recovery model and generates the HFMSE trajectory figures shown on the poster

* How to run:

1. Open a script in RStudio.
2. Install required packages listed at the top of the file.
3. Run the script to generate model predictions and figures.

* Data & Compliance:

The scripts generate model-based predictions and figures.

No patient-level or identifiable data are included.

* Reproducibility:

Running the scripts will reproduce the model outputs and visualizations presented on the poster.

* License:

MIT License

* Citation:
If you reference or build upon this work, please credit:

**Anish Kommalapati, 2026**
