import json, csv

J = json.load(open("pb1_clean.fasta.FEL.json"))
mle = J["MLE"]
heads = [h[0].lower() for h in mle["headers"]]
rows = mle["content"]["0"]
ia, ib, ip = heads.index("alpha"), heads.index("beta"), heads.index("p-value")

# PB1 C-terminus residues that contact PB2 (Sugiyama et al. 2009)
contact = {692,695,696,698,699,715,719,722,725,746,747,750,753,754}
P = 0.1
pos, neg = [], []
with open("fel_sites.csv","w",newline="") as f:
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
