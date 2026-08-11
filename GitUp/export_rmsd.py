import os, csv, math
import numpy as np
from analyze import parse_cif, analyze, FILES, DIR

OFFSET = 261  # file residue + 261 = human RGS20 numbering

# ---------------------------------------------------------------------------
# (a)+(c)  Common-contact table in HUMAN RGS20 numbering -> CSV
# ---------------------------------------------------------------------------
results = {tag: analyze(os.path.join(DIR, fn)) for tag, fn in FILES.items()}
sets = {tag: set(results[tag]["pairs"].keys()) for tag in FILES}
common = sets["rank_1"] & sets["rank_2"] & sets["rank_3"]

PR = {"Salt bridge": 0, "H-bond": 1, "Hydrophobic": 2, "van der Waals": 3}
def best_type(types):
    return sorted(types, key=lambda t: PR[t])[0]

rows = []
for key in common:
    ar, arn, br, brn = key
    ds = [results[t]["pairs"][key][0] for t in FILES]
    types = [results[t]["pairs"][key][1] for t in FILES]
    ctype = best_type(types)
    pA = results["rank_1"]["res_plddt"]["A"][ar][1]
    pB = results["rank_1"]["res_plddt"]["B"][br][1]
    pB2 = results["rank_2"]["res_plddt"]["B"][br][1]
    pB3 = results["rank_3"]["res_plddt"]["B"][br][1]
    pA2 = results["rank_2"]["res_plddt"]["A"][ar][1]
    pA3 = results["rank_3"]["res_plddt"]["A"][ar][1]
    rows.append({
        "GNAZ_residue": f"{arn}{ar}",
        "RGSZ1_residue_human": f"{brn}{br + OFFSET}",
        "RGSZ1_residue_file": f"{brn}{br}",
        "contact_type": ctype,
        "type_rank1": types[0], "type_rank2": types[1], "type_rank3": types[2],
        "consistent_type": "yes" if len(set(types)) == 1 else "no",
        "mindist_rank1": round(ds[0], 2),
        "mindist_rank2": round(ds[1], 2),
        "mindist_rank3": round(ds[2], 2),
        "mean_dist": round(sum(ds) / 3, 2),
        "pLDDT_GNAZ_mean": round((pA + pA2 + pA3) / 3, 1),
        "pLDDT_RGSZ1_mean": round((pB + pB2 + pB3) / 3, 1),
    })
rows.sort(key=lambda r: (int("".join(c for c in r["GNAZ_residue"] if c.isdigit())),
                         int("".join(c for c in r["RGSZ1_residue_file"] if c.isdigit()))))

out_csv = os.path.join(DIR, "common_contacts_RGSdomain_numbering.csv")
fields = ["GNAZ_residue", "RGSZ1_residue_human", "RGSZ1_residue_file", "contact_type",
          "type_rank1", "type_rank2", "type_rank3", "consistent_type",
          "mindist_rank1", "mindist_rank2", "mindist_rank3", "mean_dist",
          "pLDDT_GNAZ_mean", "pLDDT_RGSZ1_mean"]
with open(out_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
print(f"[a/c] wrote {len(rows)} common contacts -> {out_csv}")

# ---------------------------------------------------------------------------
# (b)  Interface Cα-RMSD between the three models (Kabsch superposition)
# ---------------------------------------------------------------------------
def ca_coords(path):
    """returns dict (chain,resnum)->np.array(xyz) for CA atoms of chains A,B"""
    atoms, _ = parse_cif(path)
    ca = {}
    for ch in ("A", "B"):
        for (rn, rname, aname, el, x, y, z, p) in atoms.get(ch, []):
            if aname == "CA":
                ca[(ch, rn)] = np.array([x, y, z])
    return ca

CA = {tag: ca_coords(os.path.join(DIR, fn)) for tag, fn in FILES.items()}

# interface residues = union of all interface residue pairs across models
iface_A = set(("A", k[0]) for t in FILES for k in sets[t])
iface_B = set(("B", k[2]) for t in FILES for k in sets[t])
iface_res = iface_A | iface_B

def kabsch_rmsd(P, Q):
    """RMSD after optimal superposition of P onto Q (Nx3 arrays)."""
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Wt.T @ V.T))
    D = np.diag([1, 1, d])
    R = Wt.T @ D @ V.T
    Prot = Pc @ R.T
    diff = Prot - Qc
    return math.sqrt((diff * diff).sum() / len(P))

def rmsd_over(keys, a, b):
    common_k = [k for k in keys if k in CA[a] and k in CA[b]]
    P = np.array([CA[a][k] for k in common_k])
    Q = np.array([CA[b][k] for k in common_k])
    return kabsch_rmsd(P, Q), len(common_k)

all_keys = set(CA["rank_1"]) | set(CA["rank_2"]) | set(CA["rank_3"])
pairs = [("rank_1", "rank_2"), ("rank_1", "rank_3"), ("rank_2", "rank_3")]

print("\n[b] Cα-RMSD between models (Kabsch best-fit superposition)")
print(f"{'pair':<16}{'all-Cα (A+B)':>16}{'interface-Cα':>16}{'RGSZ1-only':>14}{'GNAZ-only':>13}")
for a, b in pairs:
    r_all, n_all = rmsd_over(all_keys, a, b)
    r_if, n_if = rmsd_over(iface_res, a, b)
    r_b, n_b = rmsd_over(set(k for k in all_keys if k[0] == "B"), a, b)
    r_a, n_a = rmsd_over(set(k for k in all_keys if k[0] == "A"), a, b)
    print(f"{a+' vs '+b:<16}{r_all:>13.3f} A{r_if:>13.3f} A{r_b:>11.3f} A{r_a:>10.3f} A")
print(f"\n(all-Cα n={n_all}, interface-Cα n={n_if}, RGSZ1 n={n_b}, GNAZ n={n_a})")

# Domain-placement test: superpose on GNAZ only, then measure RGSZ1 displacement
print("\n[b'] RGSZ1 displacement after superposing on GNAZ core only")
def superpose_then_measure(a, b, fit_keys, meas_keys):
    fit = [k for k in fit_keys if k in CA[a] and k in CA[b]]
    P = np.array([CA[a][k] for k in fit]); Q = np.array([CA[b][k] for k in fit])
    Pc0 = P.mean(0); Qc0 = Q.mean(0)
    Pc = P - Pc0; Qc = Q - Qc0
    H = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Wt.T @ V.T))
    R = Wt.T @ np.diag([1, 1, d]) @ V.T
    meas = [k for k in meas_keys if k in CA[a] and k in CA[b]]
    Pm = (np.array([CA[a][k] for k in meas]) - Pc0) @ R.T
    Qm = np.array([CA[b][k] for k in meas]) - Qc0
    diff = Pm - Qm
    return math.sqrt((diff * diff).sum() / len(meas))
gnaz_keys = set(k for k in all_keys if k[0] == "A")
rgsz1_keys = set(k for k in all_keys if k[0] == "B")
for a, b in pairs:
    r = superpose_then_measure(a, b, gnaz_keys, rgsz1_keys)
    print(f"  {a} vs {b}: RGSZ1 Cα-RMSD after GNAZ fit = {r:.3f} A")
