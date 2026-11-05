# 2D Stability Selection

R code for the simulations and real-data analysis in
**"2D Stability Selection: Design Jittering to Exclude Correlated Irrelevant Variables"**
(Nouraie, Zhu and Muller).

2D-SS runs Stability Selection over a grid of Gaussian design-jittering levels δ
(δ = 0 is ordinary Stability Selection) and picks the selected set with a largest-gap rule.
A jittered Lasso fit is, in expectation, an elastic net, so the jittering acts as a ridge penalty.

## Folders

Each folder holds one script, named after the study, and the CSV and PDF files it writes.
The folder names match the sections of the paper.

| Folder | Script | Content in the paper |
|---|---|---|
| `Fig2` | `Fig1.R` | Motivating example (irrepresentable condition violated) |
| `S3_1` | `S31.R` | Study 1: theory curves (Corollary 2) |
| `S3_2` | `S32.R` | Study 2: selection performance and confounder exclusion |
| `S3_3` | `S33.R` | Study 3: error-bound results |
| `S3_4` | `S34.R` | Study 4: confounders, F1, sensitivity, selection rules |
| `S3_4_addon` | `S34addon.R` | Study 4: additional runs |
| `S3_5` | `S35.R` | Study 5: profiles, rules, subset size |
| `S3_6` | `S36.R` | Study 6: ENet-SS baseline (deterministic elastic net) |
| `Rat` | `Rat.R` | Rat data: repeated random splits |
| `Sensitivity` | `Sensitivity.R` | Sensitivity analysis (appendix) |
| `Diagram` | – | Figure 1 (schematic) |

Some output files keep older labels (`study0`, `corollary1`, `fig1_`).
`S3_1` is Study 1, and the Fig2 files are for the motivating example.

## Usage

1. Open the script in R and set the working directory to its folder.
2. Run the script. It writes `section3_*.csv` (results) and `plot_*.pdf` (figures) into the same folder.
3. The PDFs are the figures used in the paper.

## Contact

Mahdi Nouraie, soleil.nouraie@mq.edu.au
