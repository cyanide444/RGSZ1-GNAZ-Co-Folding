import os
from analyze import parse_cif, FILES, DIR

for tag, fn in FILES.items():
    _, plddt_res = parse_cif(os.path.join(DIR, fn))
    print(f"\n=== {tag} : residues with mean pLDDT < 70 ===")
    for ch, label in (("A","GNAZ"),("B","RGSZ1")):
        lows = []
        for rn, lst in sorted(plddt_res[ch].items()):
            m = sum(v for _,v in lst)/len(lst)
            if m < 70:
                lows.append((rn, lst[0][0], m))
        if lows:
            # collapse into ranges-ish print
            print(f"  {label}: " + ", ".join(f"{n}{r}({m:.0f})" for r,n,m in lows))
        else:
            print(f"  {label}: none < 70")
    # also report 70-80 band briefly for GNAZ termini
