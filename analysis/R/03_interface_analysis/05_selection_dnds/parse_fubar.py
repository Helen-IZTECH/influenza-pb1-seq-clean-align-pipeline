import json, csv, sys

# Usage: python3 parse_fubar.py            -> PB1 (reads pb1_clean.fasta.FUBAR.json, writes fubar_sites.csv)
#        python3 parse_fubar.py PB2 [contact_sites_csv]
#                                      -> reads pb2_clean.fasta.FUBAR.json, writes fubar_sites_PB2.csv
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

J = json.load(open(f"{PREFIX}.fasta.FUBAR.json"))
mle = J["MLE"]
heads = [h[0].lower() for h in mle["headers"]]
rows = mle["content"]["0"]
ia, ib = heads.index("alpha"), heads.index("beta")
ipur = heads.index("prob[alpha>beta]")   # purifying posterior
ipos = heads.index("prob[alpha<beta]")   # positive posterior

TH = 0.9
pos, pur = [], []
with open(f"fubar_sites{SUFFIX}.csv","w",newline="") as f:
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
