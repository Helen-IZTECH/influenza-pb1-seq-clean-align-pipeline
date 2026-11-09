import json, csv, sys

# Usage: python3 parse_fel.py            -> PB1 (reads pb1_clean.fasta.FEL.json, writes fel_sites.csv)
#        python3 parse_fel.py PB2 [contact_sites_csv]
#                                      -> reads pb2_clean.fasta.FEL.json, writes fel_sites_PB2.csv
PROTEIN = sys.argv[1].upper() if len(sys.argv) > 1 else "PB1"
PREFIX  = PROTEIN.lower() + "_clean"
SUFFIX  = "" if PROTEIN == "PB1" else "_" + PROTEIN
if PROTEIN == "PB1":
    # PB1 C-terminus residues that contact PB2 (Sugiyama et al. 2009)
    contact = {692,695,696,698,699,715,719,722,725,746,747,750,753,754}
else:
    # PB2 N-terminal residues contacting PB1 (comma list; default = 3A1G contacts, 4.5 A)
    contact = {int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else
               "1,2,3,4,5,6,7,8,10,11,15,16,17,19,20,21,23,24,25,26,28,31,32,34,35").split(",")}

J = json.load(open(f"{PREFIX}.fasta.FEL.json"))
mle = J["MLE"]
heads = [h[0].lower() for h in mle["headers"]]
rows = mle["content"]["0"]
ia, ib, ip = heads.index("alpha"), heads.index("beta"), heads.index("p-value")

P = 0.1
pos, neg = [], []
with open(f"fel_sites{SUFFIX}.csv","w",newline="") as f:
    w = csv.writer(f)
    w.writerow(["site","alpha","beta","p_value","class","is_contact"])
    for i, r in enumerate(rows, start=1):
        a, b, p = r[ia], r[ib], r[ip]
        cls = "neutral"
        if   p <= P and b > a: cls = "positive"; pos.append(i)
        elif p <= P and b < a: cls = "negative"; neg.append(i)
        w.writerow([i, round(a,4), round(b,4), round(p,4), cls, i in contact])

print("sites:", len(rows), "| positive:", pos, "| purifying:", len(neg))
print("contact under purifying:", [i for i in sorted(contact) if i in set(neg)])
print("contact under positive :", [i for i in sorted(contact) if i in set(pos)] or "none")
