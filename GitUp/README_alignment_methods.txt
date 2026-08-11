README - ALIGNMENT & CONSERVATION METHODS (1AGR restraint mapping)
==================================================================
Directory : C:\Users\richi\Documents\1AGR
Purpose   : Derive HADDOCK restraint residues for the GNAZ.RGSZ1 docking from the
            1AGR RGS4-Gai1 transition-state complex, and test whether the 1AGR
            interface residues are conserved (a) between 1AGR and the human
            orthologs and (b) in the actual docking targets GNAZ and RGSZ1.

Scripts (this directory)
  contacts_1agr.py  - interface detection + restraint assessment on 1agr_mod.pdb
  align_1agr.py     - sequence alignments, residue mapping, conservation calls
Outputs (this directory)
  _iface.json                        - interface residue lists + 1AGR chain seqs
  1AGR_alignments.txt                - the four pairwise alignments (printed)
  1AGR_conservation_GNAZ_RGSZ1.txt   - per-residue conservation of the restraints


1. INPUT STRUCTURE
------------------
1agr_mod.pdb = 1AGR (Tesmer et al., Cell 1997): RGS4 . Gai1 . GDP . AlF4- . Mg2+
transition-state complex.
  chain A = Gai1  (resolved residues 5-355; native Gai1 numbering)
  chain E = RGS4  (resolved residues 51-178; native RGS4 numbering)
  HETATM  = ALF (AlF4-), CIT, MG, HOH  (excluded from the protein-protein analysis)


2. INTERFACE DETECTION  (contacts_1agr.py)
------------------------------------------
- Only ATOM records of chains A and E are read; hydrogens (element H) are dropped.
- A Gai1 residue and an RGS4 residue are counted as an interface pair if ANY pair
  of their heavy atoms is within 5.0 A (DET = 5.0). This is the same 5.0 A
  heavy-atom criterion used throughout the project's contact analyses.
- For each interface residue pair the tightest contact is classified, in priority
  order:
    Salt bridge : cationic (Lys NZ; Arg NE/NH1/NH2) to anionic (Asp OD1/OD2;
                  Glu OE1/OE2) at <= 4.0 A
    H-bond      : both atoms are N or O and <= 3.5 A (heavy-atom donor/acceptor;
                  no explicit-H geometry, since 1AGR has no hydrogens)
    Hydrophobic : C-C contact <= 4.5 A between two hydrophobic residues
                  (A,V,L,I,M,F,W,P,C)
    van der Waals: any remaining contact <= 5.0 A
  A residue's reported type is the strongest type over all its cross-interface
  contacts.
- ACTIVE vs PASSIVE split: an interface residue whose strongest contact is a salt
  bridge or H-bond is called active; a residue that only makes hydrophobic/vdW
  contacts is called passive. (Active = directly-interacting anchor; passive =
  surrounding surface, HADDOCK convention.)
- The user-proposed residue lists are compared to the computed interface to flag
  (i) proposed residues that do NOT contact (candidates for removal) and
  (ii) interface residues absent from the proposal (candidates to add).
- Interface residue lists and the extracted 1AGR chain sequences are written to
  _iface.json for the alignment step.


3. SEQUENCE ALIGNMENT  (align_1agr.py)
--------------------------------------
Algorithm: global Needleman-Wunsch, BLOSUM62 substitution matrix, linear gap
penalty. Full dynamic-programming matrix with traceback; no heuristics. Match
markers in the printed alignments: '|' identical, ':' non-identical but positive
BLOSUM62 score, ' ' otherwise.

Reference sequences (UniProt canonical, full length):
  human Gai1  P63096
  human RGS4  P49798
  GNAZ (Gai-z) P19086
  RGSZ1 RGS domain = the docked fragment (residues 1-121 as numbered in
        RGSZ1_protein_only_docking.pdb; equals human RGSZ1/RGS20 minus 259)
1AGR chain sequences are reconstructed from the resolved residues in _iface.json.

Four alignments are produced (1AGR_alignments.txt):
  1. 1AGR Gai1 (chain A, resolved)  vs  human Gai1 (P63096)
  2. human Gai1 (P63096)            vs  GNAZ (P19086)
  3. 1AGR RGS4 (chain E, resolved)  vs  human RGS4 (P49798)
  4. human RGS4 (P49798)            vs  RGSZ1 RGS domain (docked 1-121)

Gap-penalty note for alignment 4:
  Human RGS4 carries a ~50-residue N-terminal extension that has no counterpart in
  the RGSZ1 1-121 domain fragment. A default global alignment spreads spurious gaps
  through the termini. Alignment 4 therefore aligns only the RGS-box core (human
  RGS4 from residue 51) with a stiffer gap penalty (-11 vs the default -8) so the
  conserved core aligns cleanly. Even so, the RGS4 alpha6-alpha7 loop (around
  R166/R167) genuinely differs in length between RGS4 and RGSZ1 (a real indel), so
  those two positions receive no 1:1 map - this is a biological divergence, not an
  alignment artifact.


4. RESIDUE MAPPING & CONSERVATION CALLS
---------------------------------------
- From each pairwise alignment a residue-number map is built: for every aligned
  (non-gap, non-gap) column, source residue number -> (target residue, target
  number). Residue numbers are anchored to the true start of each sequence
  (alignment 4 uses base 51 for RGS4).
- 1AGR vs human: because 1AGR uses native Gai1 / RGS4 numbering, conservation
  between 1AGR and human is read by direct same-number comparison; every proposed
  and interface residue was identical (1AGR = human).
- Target mapping: human Gai1 -> GNAZ gives a uniform +1 offset (Gai1 N -> GNAZ N+1);
  human RGS4 -> RGSZ1 gives a -59 core offset (RGS4 N -> RGSZ1 N-59).
- Conservation label for a restraint residue in the docking target:
    conserved = target residue identity equals the human/1AGR residue identity
    similar   = mapped but different residue (substitution flagged)
    no map    = falls in a gap / divergent loop (no 1:1 residue)
  Charge reversals and size changes are called out explicitly (e.g. Gai1 E116 ->
  GNAZ K117; RGS4 N128 -> RGSZ1 S69).


5. LIMITATIONS
--------------
- 1AGR contains no hydrogens, so H-bonds are assigned by heavy-atom distance
  (<= 3.5 A) without donor-H-acceptor angle checks.
- Alignments are pairwise/global, not a curated structural superposition; the
  alpha6-alpha7 RGS4/RGSZ1 indel is reported as unmapped rather than force-fit.
- Conservation is sequence-identity based; a non-identical but functionally
  equivalent substitution is labeled "similar", not "conserved".


REPRODUCE
---------
  cd C:\Users\richi\Documents\1AGR
  py -3 contacts_1agr.py     # writes _iface.json
  py -3 align_1agr.py        # writes 1AGR_alignments.txt and
                             #        1AGR_conservation_GNAZ_RGSZ1.txt
