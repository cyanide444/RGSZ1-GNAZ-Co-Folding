# per_residue_deviation.py
#
# Called from condense_ranks.cxc after the three AF3 rank models have been
# opened and superimposed. For each residue (matched by chain ID + auth
# residue number, which is consistent across rank_1/2/3.cif), this computes
# the RMS spread of the CA atom position across the three ranks and stores
# it as a custom residue attribute "rankdev" on model #1. That model can
# then be colored/painted by this attribute to show where the three ranks
# disagree.
#
# Requires: the three models already open in the session as #1, #2, #3
# (adjust the id strings below if your models load with different numbers)
# and already superimposed (matchmaker in the .cxc handles this).

import numpy as np
from chimerax.atomic import AtomicStructure, Residue

# Register the custom attribute so `color byattribute` / `select` can see it
Residue.register_attr(session, "rankdev", "rgsz1_gnaz_script", attr_type=float)

mdls = {m.id_string: m for m in session.models if isinstance(m, AtomicStructure)}
missing = [i for i in ("1", "2", "3") if i not in mdls]
if missing:
    print(f"per_residue_deviation.py: could not find model(s) {missing} — "
          f"check the model numbers in ChimeraX (see the Log/Models panel) "
          f"and edit the id strings in this script if needed.")
else:
    m1, m2, m3 = mdls["1"], mdls["2"], mdls["3"]

    def get_ca(m, chain_id, resnum):
        for res in m.residues:
            if res.chain_id == chain_id and res.number == resnum:
                a = res.find_atom("CA")
                if a is not None:
                    return a.scene_coord, res
        return None, None

    # chain A = GNAZ (1-355), chain B = RGSZ1 fragment (1-117, = true 262-378)
    chain_ranges = {"A": range(1, 356), "B": range(1, 118)}

    n_scored = 0
    max_dev = 0.0
    for chain_id, rng in chain_ranges.items():
        for resnum in rng:
            p1, r1 = get_ca(m1, chain_id, resnum)
            p2, _ = get_ca(m2, chain_id, resnum)
            p3, _ = get_ca(m3, chain_id, resnum)
            pts = [p for p in (p1, p2, p3) if p is not None]
            if len(pts) < 2 or r1 is None:
                continue
            pts = np.array(pts)
            centroid = pts.mean(axis=0)
            dev = float(np.sqrt(((pts - centroid) ** 2).sum(axis=1).mean()))
            r1.rankdev = dev
            n_scored += 1
            max_dev = max(max_dev, dev)

    print(f"per_residue_deviation.py: scored {n_scored} residues on model #1. "
          f"Max per-residue deviation = {max_dev:.2f} \u00c5. "
          f"Now run: color byattribute r:rankdev #1 palette blue:white:red")
