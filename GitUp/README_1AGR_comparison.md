# Comparison to the RGS4·Gαi1 transition-state complex (1AGR)

Places the co-folded GNAZ·RGSZ1 model in the context of the only high-resolution RGS·Gα
complex, **1AGR** (RGS4·Gαi1·GDP·AlF₄⁻·Mg²⁺, a transition-state mimic; rat sequences,
chain A = Gαi1, chain E = RGS4). Establishes the residue correspondence, compares the two
interfaces, and evaluates how transition-state-like the modeled active site is.

See also `README_alignment_methods.txt` (original method notes) and
`RGSZ1_RGS_domain_topology.txt` (RGSZ1 → RGS-box helix map).

## Scripts

| Script | Role | Output |
|--------|------|--------|
| `align_1agr.py` | Pairwise sequence alignment (BLOSUM62; Needleman–Wunsch global and Smith–Waterman local) of **GNAZ vs Gαi1** and **RGSZ1 vs RGS4**, to map homologous residues and score conservation across the interface. | `1AGR_alignments.txt`, `1AGR_conservation_GNAZ_RGSZ1.txt` |
| `contacts_1agr.py` | Computes the RGS4–Gαi1 interface contacts in 1AGR (same geometric definitions as the co-fold contact analysis) so the two interfaces can be compared contact-for-contact. | interface JSON / conservation table |
| `map_helices.py` | Maps RGSZ1 residues onto the nine Tesmer et al. RGS-box helices via the RGS4 homology (RGS4 residue = RGSZ1 file residue + 61). | `RGSZ1_RGS_domain_topology.txt` |
| `ts_analysis.py` | Reads the ChimeraX Matchmaker superposition of **1AGR vs rank_1** (`Structures/_1agr_super.pdb`, `Structures/_rank1_super.pdb`) and measures the active-site geometry: catalytic Gln (Gαi1 Gln204 / GNAZ Gln205), the arginine finger, the γ-phosphate/Mg²⁺, and the catalytic Asn128(RGS4)/Ser328(RGSZ1) — reporting how closely the modeled site reproduces the transition state. | `TS_analysis_rank1_vs_1AGR.txt` |

## Key points captured by these analyses

- **Residue correspondence.** The RGS4 catalytic **Asn128** aligns to RGSZ1 **Ser328**
  (file Ser67) — the RZ-subfamily Asn→Ser substitution. In 1AGR, Asn128 hydrogen-bonds the
  catalytic Gln204 of Gαi1 (~3.0 Å); in the model the RGSZ1 serine does not reach GNAZ Gln205
  and instead engages Switch II Glu208.
- **Transition-state character.** After superposition, the GNAZ Gln205 side chain, the
  γ-phosphate/Mg²⁺ and the arginine finger overlay their 1AGR counterparts within ~0.3–0.9 Å,
  i.e. the modeled active site is pre-organized in a transition-state-like geometry
  (`TS_analysis_rank1_vs_1AGR.txt`).
- **Interface comparison.** `contacts_1agr.py` + `align_1agr.py` show which RGS4·Gαi1 contacts
  are conserved, substituted, or absent in the GNAZ·RGSZ1 model.

## Structures used
`Structures/1AGR.pdb`, and the superimposed coordinates `Structures/_rank1_super.pdb` /
`Structures/_1agr_super.pdb`. The ChimeraX sessions `rank1(4contact)_1AGRsuperimposed.cxs`
and `1agr_gln204-asn128.cxs` open the superposition and the active-site view directly.
