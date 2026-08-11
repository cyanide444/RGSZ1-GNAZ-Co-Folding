# Interface contact analysis (GNAZ · RGSZ1)

Identifies and classifies the residue–residue contacts across the GNAZ(chain A)–RGSZ1(chain B)
interface in the three co-folded models, and reports the set that is **reproducible across all
three ranks** as the confident interface.

## Method

**Contact definition.** Heavy atoms only (hydrogens ignored except where ChimeraX adds them
for H-bond geometry). Two residues are in contact if any heavy-atom pair is within the cutoff.

**Contact typing** (a residue pair is assigned the single strongest type present), in priority
order:
1. **Salt bridge** — cationic N (Arg NH1/NH2/NE, Lys NZ, His ND1/NE2) ↔ anionic O
   (Asp OD1/OD2, Glu OE1/OE2) within **4.0 Å**.
2. **Hydrogen bond** — taken from **ChimeraX** (`addh` then `hbonds`), i.e. proper
   donor/acceptor identity + distance + angle geometry, using the relaxed default tolerances
   (0.4 Å / 20°). Precomputed into `hb_rank*.txt` / `hb_cif*.txt` by `chimerax_hbonds.cxc`.
3. **Hydrophobic** — C–C contact ≤ **4.5 Å** between two hydrophobic residues.
4. **van der Waals** — any heavy-atom pair ≤ **4.0 Å** not already classified.

**Reproducibility filter.** A contact is "common" if present in **rank_1 AND rank_2 AND rank_3**.
The `_filtered` table additionally requires both partner residues to have mean **pLDDT ≥ 85**
(minimum across the three models), i.e. the contact lies in a confidently modeled region.

**Region labels.** GNAZ residues are labeled Switch I / II / III / other; RGSZ1 residues are
labeled by Tesmer et al. RGS-box helix/loop via the RGS4 homology (RGS4 residue = RGSZ1 file
residue + 61). RGSZ1 residues are reported in both file (1–117) and human (262–378) numbering.

## Scripts

| Script | Role |
|--------|------|
| `chimerax_hbonds.cxc` | ChimeraX: `open` rank_1/2/3, `addh`, `hbonds #n/A restrict #n/B` → `hb_rank{1,2,3}.txt`. Run this first. |
| `analyze.py` | Core module. Parses `rank_*.cif`, computes typed contacts and per-residue pLDDT. Imported by the two scripts below. |
| `contacts.py` | Standalone quick look — lists heavy-atom interface contacts (≤ 4.0 Å) for a single model and prints a chain/residue summary. |
| `reassess.py` | Rebuilds the common-contact set from the **PDB** files, merging the ChimeraX H-bonds with the geometric salt-bridge/hydrophobic/vdW calls; applies the pLDDT ≥ 85 filter → **`common_contacts_filtered.csv`**. |
| `regen_cif.py` | Rebuilds the table directly from the **CIF** files (authoritative `auth_seq_id` numbering) and adds the Switch / Tesmer-helix region labels. |
| `export_rmsd.py` | Imports `analyze`; writes the common-contact table in human RGS20 numbering with per-rank distances, per-rank types, consistency flag and pLDDT. |
| `lowconf.py` | Imports `analyze`; lists residues with mean pLDDT < 70 per chain per model (confidence QC). |
| `per_residue_deviation.py` | ChimeraX script: after the three ranks are superimposed, computes the per-residue Cα spread across ranks and stores it as the residue attribute `rankdev` for coloring. |
| `condense_ranks.cxc` | ChimeraX driver: opens + Matchmaker-superimposes the three ranks and runs `per_residue_deviation.py` (source of the Cα-deviation figures). |

## Outputs (included)

- **`common_contacts_filtered.csv`** — the confident, reproducible interface (27 contacts;
  columns: GNAZ residue + region, RGSZ1 residue in human and file numbering, Tesmer element,
  contact type, per-rank min distances, mean distance, consistency).
- **`common_contacts_RGSdomain_numbering.csv`** — same contacts keyed on RGS-domain numbering.

The strongest single contact is the **Lys211(GNAZ)–Glu287(RGSZ1)** salt bridge; the catalytic
region contact **Glu208(GNAZ)–Ser328(RGSZ1)** is the one whose classification varies between
ranks (it is the RZ-specific serine described in the 1AGR comparison).
