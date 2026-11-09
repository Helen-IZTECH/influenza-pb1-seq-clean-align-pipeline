import json, csv, sys

# Usage: python3 parse_meme.py            -> PB1 (reads pb1_clean.fasta.MEME.json, writes meme_sites.csv)
#        python3 parse_meme.py PB2 [contact_sites_csv]
#                                      -> reads pb2_clean.fasta.MEME.json, writes meme_sites_PB2.csv
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

J = json.load(open(f"{PREFIX}.fasta.MEME.json"))
mle = J["MLE"]
heads = [h[0].lower() for h in mle["headers"]]
rows = mle["content"]["0"]
ip = [i for i, h in enumerate(heads) if "p-value" in h][0]

P = 0.1
sig = []
with open(f"meme_sites{SUFFIX}.csv","w",newline="") as f:
    w = csv.writer(f)
    w.writerow(["site","p_value","episodic_positive","is_contact"])
    for i, r in enumerate(rows, start=1):
        p = r[ip]
        e = p <= P
        if e: sig.append(i)
        w.writerow([i, round(p,4), e, i in contact])

print("episodic positive (p<=0.1):", sig)
print("contact among them:", [i for i in sorted(contact) if i in set(sig)] or "none")
