import json, csv

J = json.load(open("pb1_clean.fasta.FUBAR.json"))
mle = J["MLE"]
heads = [h[0].lower() for h in mle["headers"]]
rows = mle["content"]["0"]
ia, ib = heads.index("alpha"), heads.index("beta")
ipur = heads.index("prob[alpha>beta]")   # purifying posterior
ipos = heads.index("prob[alpha<beta]")   # positive posterior

contact = {692,695,696,698,699,715,719,722,725,746,747,750,753,754}
TH = 0.9
pos, pur = [], []
with open("fubar_sites.csv","w",newline="") as f:
    w = csv.writer(f)
    w.writerow(["site","alpha","beta","purifying_post","positive_post","is_contact"])
    for i, r in enumerate(rows, start=1):
        a, b, pp, po = r[ia], r[ib], r[ipur], r[ipos]
        if po >= TH: pos.append(i)
        if pp >= TH: pur.append(i)
        w.writerow([i, round(a,3), round(b,3), round(pp,3), round(po,3), i in contact])

print("sites:", len(rows), "| positive>=0.9:", pos, "| purifying>=0.9:", len(pur))
print("contact positive :", [i for i in sorted(contact) if i in set(pos)] or "none")
print("contact purifying:", [i for i in sorted(contact) if i in set(pur)])
