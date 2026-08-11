import math, os

DIR = r"C:\Users\richi\Documents\RGSZ1 Human (262-378)-GNAZ-GTP Complex 4contact"
FILES = {"rank_1": "rank_1.cif", "rank_2": "rank_2.cif", "rank_3": "rank_3.cif"}
CUTOFF = 4.0          # heavy-atom contact cutoff (vdW/general)
HBOND = 3.5           # polar N/O - N/O
SALT = 4.0            # charged - charged
HPHOB = 4.5           # C - C apolar

# charged atom sets
POS = {("ARG", a) for a in ("NH1", "NH2", "NE")} | {("LYS", "NZ")} | {("HIS", "ND1"), ("HIS", "NE2")}
NEG = {("ASP", a) for a in ("OD1", "OD2")} | {("GLU", a) for a in ("OE1", "OE2")}
DONOR_ACC = set("NO")  # element-based polar
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS", "TYR"}

def parse_cif(path):
    """returns dict chain -> list of atoms (resnum,resname,atomname,elem,x,y,z,plddt)"""
    atoms = {}
    plddt_res = {}  # chain -> {resnum: [plddt values per atom]}
    with open(path) as fh:
        in_loop = False
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                f = line.split()
                if len(f) < 19:
                    continue
                rec = f[0]
                elem = f[2]
                atomname = f[3]
                resname = f[5]
                chain = f[15]          # auth_asym_id
                try:
                    resnum = int(f[7])  # auth_seq_id
                    x = float(f[10]); y = float(f[11]); z = float(f[12])
                    plddt = float(f[17])
                except ValueError:
                    continue
                if elem == "H":
                    continue
                atoms.setdefault(chain, []).append(
                    (resnum, resname, atomname, elem, x, y, z, plddt))
                plddt_res.setdefault(chain, {}).setdefault(resnum, []).append((resname, plddt))
    return atoms, plddt_res

def classify(an_res, an_atom, an_el, bn_res, bn_atom, bn_el, d):
    """classify a single closest-atom contact between residue pair"""
    a_key = (an_res, an_atom); b_key = (bn_res, bn_atom)
    # salt bridge
    if d <= SALT and ((a_key in POS and b_key in NEG) or (a_key in NEG and b_key in POS)):
        return "Salt bridge"
    # hydrogen bond (polar N/O - N/O within hbond dist)
    if d <= HBOND and an_el in DONOR_ACC and bn_el in DONOR_ACC:
        return "H-bond"
    # hydrophobic C-C
    if d <= HPHOB and an_el == "C" and bn_el == "C" \
            and an_res in HYDROPHOBIC_RES and bn_res in HYDROPHOBIC_RES:
        return "Hydrophobic"
    if d <= CUTOFF:
        return "van der Waals"
    return None

def analyze(path):
    atoms, plddt_res = parse_cif(path)
    A = atoms.get("A", []); B = atoms.get("B", [])
    GTP = atoms.get("C", []); MG = atoms.get("D", [])

    # residue-level mean pLDDT
    res_plddt = {}
    for ch in ("A", "B"):
        res_plddt[ch] = {}
        for rn, lst in plddt_res.get(ch, {}).items():
            res_plddt[ch][rn] = (lst[0][0], sum(v for _, v in lst) / len(lst))

    # A-B interface contacts: keep per residue-pair the best (closest) contact + best classification
    cut2 = CUTOFF * CUTOFF
    # store all qualifying atom-atom for classification priority
    pair_best = {}   # (ar,arn,br,brn) -> (mindist, type, aatom, batom)
    PRIORITY = {"Salt bridge": 0, "H-bond": 1, "Hydrophobic": 2, "van der Waals": 3}
    for (ar, arn, aa, ael, ax, ay, az, ap) in A:
        for (br, brn, ba, bel, bx, by, bz, bp) in B:
            dx = ax-bx; dy = ay-by; dz = az-bz
            d2 = dx*dx+dy*dy+dz*dz
            if d2 <= HPHOB*HPHOB:
                d = math.sqrt(d2)
                t = classify(arn, aa, ael, brn, ba, bel, d)
                if t is None:
                    continue
                key = (ar, arn, br, brn)
                cur = pair_best.get(key)
                # prefer stronger type; within same, prefer shorter dist
                if cur is None or PRIORITY[t] < PRIORITY[cur[1]] or \
                   (PRIORITY[t] == PRIORITY[cur[1]] and d < cur[0]):
                    pair_best[key] = (d, t, aa, ba)

    # ligand contacts (GTP / Mg to protein) within 4.0
    def lig_contacts(lig, name):
        res = {}
        for (lr, lrn, la, lel, lx, ly, lz, lp) in lig:
            for ch, atomlist in (("A", A), ("B", B)):
                for (ar, arn, aa, ael, ax, ay, az, ap) in atomlist:
                    d2 = (lx-ax)**2+(ly-ay)**2+(lz-az)**2
                    if d2 <= cut2:
                        d = math.sqrt(d2)
                        key = (ch, ar, arn)
                        if key not in res or d < res[key][0]:
                            res[key] = (d, la, aa)
        return res

    gtp_c = lig_contacts(GTP, "GTP")
    mg_c = lig_contacts(MG, "MG")

    return {
        "pairs": pair_best, "res_plddt": res_plddt,
        "gtp": gtp_c, "mg": mg_c,
        "natoms": {k: len(v) for k, v in atoms.items()},
    }

results = {}
for tag, fn in FILES.items():
    results[tag] = analyze(os.path.join(DIR, fn))

# ---------- per-model summary ----------
print("="*78)
print("GLOBAL pLDDT SUMMARY (mean per chain)")
print("="*78)
for tag in FILES:
    rp = results[tag]["res_plddt"]
    for ch, label in (("A", "GNAZ"), ("B", "RGSZ1")):
        vals = [v for _, v in rp[ch].values()]
        print(f"{tag}  chain {ch} ({label:6s}): mean pLDDT = {sum(vals)/len(vals):5.1f}  "
              f"min = {min(vals):5.1f}  max = {max(vals):5.1f}  (n={len(vals)} res)")
    print()

# ---------- common A-B contacts ----------
print("="*78)
print("COMMON A-B (GNAZ-RGSZ1) RESIDUE CONTACTS ACROSS ALL THREE MODELS")
print("="*78)
sets = {tag: set(results[tag]["pairs"].keys()) for tag in FILES}
common = sets["rank_1"] & sets["rank_2"] & sets["rank_3"]
print(f"rank_1 pairs={len(sets['rank_1'])}  rank_2={len(sets['rank_2'])}  "
      f"rank_3={len(sets['rank_3'])}  COMMON(all3)={len(common)}\n")

hdr = f"{'GNAZ (A)':<14}{'RGSZ1 (B)':<14}{'type (consensus)':<16}" \
      f"{'d1':>6}{'d2':>6}{'d3':>6}  {'plddt_A/plddt_B (r1)':<22}"
print(hdr); print("-"*len(hdr))
def best_type(types):
    pr = {"Salt bridge":0,"H-bond":1,"Hydrophobic":2,"van der Waals":3}
    return sorted(types, key=lambda t: pr[t])[0]
rows = []
for key in common:
    ar, arn, br, brn = key
    ds = [results[t]["pairs"][key][0] for t in FILES]
    types = [results[t]["pairs"][key][1] for t in FILES]
    ctype = best_type(types)
    pA = results["rank_1"]["res_plddt"]["A"][ar][1]
    pB = results["rank_1"]["res_plddt"]["B"][br][1]
    rows.append((ar, br, arn, brn, ctype, ds, types, pA, pB))
rows.sort(key=lambda r: (r[0], r[1]))
for ar, br, arn, brn, ctype, ds, types, pA, pB in rows:
    consensus = "" if len(set(types))==1 else " *varies:"+ "/".join(types)
    print(f"{arn+str(ar):<14}{brn+str(br):<14}{ctype:<16}"
          f"{ds[0]:6.2f}{ds[1]:6.2f}{ds[2]:6.2f}  {pA:5.1f}/{pB:5.1f}{consensus}")

# ---------- contacts unique to subsets (deviations) ----------
print("\n" + "="*78)
print("MODEL-SPECIFIC / VARIABLE CONTACTS (present in some but not all models)")
print("="*78)
allkeys = sets["rank_1"] | sets["rank_2"] | sets["rank_3"]
for key in sorted(allkeys - common):
    ar, arn, br, brn = key
    present = [t for t in FILES if key in sets[t]]
    info = []
    for t in FILES:
        if key in sets[t]:
            d, ty, _, _ = results[t]["pairs"][key]
            info.append(f"{t}={ty}({d:.2f})")
        else:
            info.append(f"{t}=-")
    print(f"{arn+str(ar):<12}{brn+str(br):<12}  " + "  ".join(info))

# ---------- interface residue lists ----------
print("\n" + "="*78)
print("INTERFACE RESIDUES PER MODEL")
print("="*78)
for tag in FILES:
    ga = sorted(set((k[0],k[1]) for k in sets[tag]))
    gb = sorted(set((k[2],k[3]) for k in sets[tag]))
    print(f"\n{tag}:")
    print("  GNAZ : " + ", ".join(f"{n}{r}" for r,n in ga))
    print("  RGSZ1: " + ", ".join(f"{n}{r}" for r,n in gb))

# ---------- GTP / Mg contacts (consensus) ----------
print("\n" + "="*78)
print("GTP (LIG1) CONTACTS - residues within 4.0 A in ALL three models")
print("="*78)
gsets = {t: set(results[t]["gtp"].keys()) for t in FILES}
gcommon = gsets["rank_1"] & gsets["rank_2"] & gsets["rank_3"]
for key in sorted(gcommon):
    ch, r, rn = key
    ds = [results[t]["gtp"][key][0] for t in FILES]
    lab = "GNAZ" if ch=="A" else "RGSZ1"
    print(f"  {lab:6s} {rn}{r:<6}  d(r1/r2/r3)={ds[0]:.2f}/{ds[1]:.2f}/{ds[2]:.2f}")

print("\nMg (LIG2) coordinating residues within 4.0 A:")
for t in FILES:
    keys = sorted(results[t]["mg"].keys())
    s = ", ".join(f"{('GNAZ' if c=='A' else 'RGSZ1')}:{rn}{r}({results[t]['mg'][(c,r,rn)][0]:.2f})" for c,r,rn in keys)
    print(f"  {t}: {s if s else '(none within 4.0 A)'}")
