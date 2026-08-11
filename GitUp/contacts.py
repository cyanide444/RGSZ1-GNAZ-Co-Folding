import math, sys

PDB = r"C:\Users\richi\Documents\RGSZ1 Human (262-378)-GNAZ-GTP Complex 4contact\rank_1.pdb"
CUTOFF = 4.0  # heavy-atom contact distance in Angstrom

def parse(path):
    atoms = {}  # chain -> list of (resnum, resname, atomname, element, x,y,z)
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            chain = line[21]
            resname = line[17:20].strip()
            # skip ligands/ions/water for protein-protein analysis
            if rec == "HETATM" and resname in ("GTP", "GDP", "MG", "HOH", "NA", "CL", "ZN"):
                continue
            atomname = line[12:16].strip()
            element = line[76:78].strip() or atomname[0]
            if element == "H":
                continue  # heavy atoms only
            try:
                resnum = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            atoms.setdefault(chain, []).append((resnum, resname, atomname, x, y, z))
    return atoms

atoms = parse(PDB)
chains = sorted(atoms.keys())
print("Chains found:", chains)
for c in chains:
    rs = sorted(set((a[0], a[1]) for a in atoms[c]))
    print("  Chain %s: %d atoms, %d residues (%d-%d)" % (
        c, len(atoms[c]), len(rs), rs[0][0], rs[-1][0]))

A = atoms.get("A", [])
B = atoms.get("B", [])

cut2 = CUTOFF * CUTOFF
pairs = {}  # (Aresnum,Aresname,Bresnum,Bresname) -> min dist
for (ar, an, aa, ax, ay, az) in A:
    for (br, bn, ba, bx, by, bz) in B:
        d2 = (ax-bx)**2 + (ay-by)**2 + (az-bz)**2
        if d2 <= cut2:
            d = math.sqrt(d2)
            key = (ar, an, br, bn)
            if key not in pairs or d < pairs[key][0]:
                pairs[key] = (d, aa, ba)

print("\n=== Inter-chain residue contacts (heavy atoms <= %.1f A) ===" % CUTOFF)
print("%-18s %-18s %8s  %s" % ("GNAZ (chain A)", "RGS20 (chain B)", "minDist", "closest atoms"))
for key in sorted(pairs.keys()):
    ar, an, br, bn = key
    d, aa, ba = pairs[key]
    print("%-4s %-12s   %-4s %-12s   %6.2f A   %s--%s" % (
        an, ar, bn, br, d, aa, ba))
print("\nTotal contacting residue pairs:", len(pairs))

ga = sorted(set((k[0], k[1]) for k in pairs))
gb = sorted(set((k[2], k[3]) for k in pairs))
print("\nGNAZ interface residues (%d):" % len(ga), ", ".join("%s%d" % (n, r) for r, n in ga))
print("RGS20 interface residues (%d):" % len(gb), ", ".join("%s%d" % (n, r) for r, n in gb))
