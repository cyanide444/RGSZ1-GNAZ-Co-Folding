# AlphaFold3 / Boltz-2 confidence plots

Reconstructs the model-confidence figures for the co-folded complex from the raw
`Inputs/scores.json` output.

## Script
`plot_scores.py` — reads `scores.json` and renders four figures (one per confidence metric,
all three ranks):

| Plot | What it shows |
|------|---------------|
| `pLDDT_plot.png` | Per-token pLDDT along the sequence, with chain dividers (GNAZ \| RGSZ1 \| GTP \| Mg) and the standard confidence bands (very high ≥90, confident 70–90, low 50–70, very low <50). |
| `PAE_plot.png` | Predicted Aligned Error matrix (per model) — inter-domain / inter-chain positional confidence. |
| `PDE_plot.png` | Predicted Distance Error matrix (per model). |
| `ipTM_plot.png` | Per-chain (pairwise) ipTM matrix and the overall ipTM per model — the interface-confidence metric. |

Token layout used for the axis dividers: **GNAZ 1–355 | RGSZ1 356–472 | GTP 473–504 | Mg 505**
(505 tokens total).

## To run
Point the input path in `plot_scores.py` at `Inputs/scores.json` (the script as-run read it
from a Downloads folder). Requires `numpy` + `matplotlib`. Outputs are written next to the
script; the rendered copies are in `Figures/`.
