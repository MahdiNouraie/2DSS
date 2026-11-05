# 2D Stability Selection

R code for the simulations and real-data analysis in
**"2D Stability Selection: Design Jittering to Exclude Correlated Irrelevant Variables"**
(Nouraie, Zhu and Muller).

2D-SS runs Stability Selection over a grid of Gaussian design-jittering levels δ
(δ = 0 is ordinary Stability Selection) and picks the selected set with a largest-gap rule.
A jittered Lasso fit is, in expectation, an elastic net, so the jittering acts as a ridge penalty.

## Folders

Each folder holds one script and the CSV and PDF files it writes. The folder names match the numerical studies in Section 3 of the paper.

| Folder | Script | Content in the paper |
|---|---|---|
| `Fig2` | `Fig1.R` | Motivating example (Section 1.3) |
| `S3_1` | `S31.R` | Study 1: Corollary 2 in a single fit |
| `S3_2` | `S32.R` | Study 2: robustness of the aggregated procedure |
| `S3_3` | `S33.R` | Study 3: does the PFER bound hold? |
| `S3_4` | `S34.R` | Study 4: comparison with measurement-error baselines |
| `S3_4_addon` | `S34addon.R` | Study 4: role of the selection rule (largest-gap, known s, fixed τ) |
| `S3_5` | `S35.R` | Study 5: when the largest-gap rule fails |
| `S3_6` | `S36.R` | Study 6: is the random noise needed? (ENet-SS) |
| `Rat` | `Rat.R` | Rat data: full-data selection and repeated train–test splits |
| `Sensitivity` | `Sensitivity.R` | Appendix: sensitivity to the noise grid and the weights |

## Usage

1. Open the script in R and set the working directory to its folder.
2. Run the script. It writes `section3_*.csv` (results) and `plot_*.pdf` (figures) into the same folder.
3. The PDFs are the figures used in the paper.

Notes:
- `Rat.R` downloads the data from GEO (series GSE5680), so it needs an internet connection.
- `S34addon.R` uses `section3_study4_summary.csv` from `S3_4` if that file is in the working directory. Otherwise it uses the values stored in the script.

## Contact

Mahdi Nouraie, soheil.nouraie@mq.edu.au
