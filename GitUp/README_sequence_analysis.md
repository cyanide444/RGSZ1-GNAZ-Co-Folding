# Sequence analyses (Gα family, RGS family) and ligand prep

Supporting sequence comparisons behind the supplementary figures, plus the GTP ligand
preparation used for co-folding. All sequences are UniProt full-length FASTA
(`Inputs/Sequences/`, retrieved via `curl https://rest.uniprot.org/uniprotkb/<ACC>.fasta`).
Alignments use an inline BLOSUM62 Needleman–Wunsch (global) implementation; percent identity
is computed over aligned (ungapped) columns.

## Scripts

| Script | Produces | Description |
|--------|----------|-------------|
| `seqfig.py` | `supplementary_figure1.png` | Primary sequences of human GNAZ, RGSZ1 and Gαi1, annotated with the GNAZ switch regions, the RGSZ1 RGS box and Tesmer helices, and the Gα α-helical-domain shading. |
| `ident.py` | identity matrices (used by `heatmaps.py`) | All-vs-all pairwise % identity for (a) five human Gα subunits — GNAZ P19086, GNAI1 P63096, GNAI2 P04899, GNAI3 P08754, GNAO P09471; and (b) the RGS boxes of all 17 human R4/R7/RZ RGS proteins. Also prints the RGS20 (RGSZ1)-vs-all ranking. |
| `heatmaps.py` | `galpha_identity_heatmap.png`, `rgs_identity_heatmap.png` | Renders the two % identity heatmaps (Gα 5×5; RGS 17×17 grouped by subfamily). |
| `rgsfig.py` | `supplementary_figure3.png` | Alignment of all 17 human R4/R7/RZ RGS boxes projected onto the RGS4 box; the eight RGS-fold helices are colored and the catalytic column (RGS4 Asn128) is boxed in red — showing Asn is invariant in R4/R7 but a Ser in all three RZ members. |
| `gz_gi_rat_fig.py` | `GNAZ_GNAI1_alignment.png` | Alignment of human GNAZ with human and rat Gαi1 (P19086 / P63096 / P10824), with the three switch regions shaded. Rat Gαi1 is the 1AGR sequence. |
| `v3k_to_v2k.py` | `Inputs/ionizedGTP_2.sdf` | Converts the ionized-GTP ligand SDF from V3000 to V2000 format for co-folding; per-atom naming recorded in `Inputs/ionizedGTP_atommap.txt`. |

## UniProt accessions (in `Inputs/Sequences/`)

- **Gα:** P19086 (GNAZ), P63096 (GNAI1 human), P10824 (GNAI1 rat), P04899 (GNAI2),
  P08754 (GNAI3), P09471 (GNAO).
- **RGS R4:** Q08116, P41220, P49796, P49798, O15539, P57771, O14921, O15492, Q9NS28, Q2M5E4.
- **RGS R7:** P49758, P49802, O75916, O94810.
- **RGS RZ:** Q9UGC6 (RGS17), P49795 (RGS19), O76081 (RGS20 / RGSZ1).

## Headline results
- GNAZ is the most divergent Gi/o α subunit — 67% identical to Gαi1/Gαi3, 60% to Gαo.
- The three RZ RGS domains form a tight block (mean 74% identity); RGS20 is closest to
  RGS19 (77%) and RGS17 (76%) and only ~36–50% identical to any R4/R7 domain.
- The catalytic Asn (RGS4 Asn128) is invariant across R4 (10/10) and R7 (4/4) but is a serine
  in all three RZ members — the substitution central to this project.

> **Note on paths:** `rgsfig.py` and `ident.py`/`heatmaps.py`/`gz_gi_rat_fig.py` read their
> FASTA inputs from the directory the script sits in (`rgsfam` / `galpha` in the original run).
> To re-run from this repo, set the `D=`/`GA=` directory constant to `Inputs/Sequences/`.
