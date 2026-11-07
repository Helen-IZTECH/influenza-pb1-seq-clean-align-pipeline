import json, csv

J = json.load(open("pb1_clean.fasta.MEME.json"))
mle = J["MLE"]
heads = [h[0].lower() for h in mle["headers"]]
rows = mle["content"]["0"]
ip = [i for i, h in enumerate(heads) if "p-value" in h][0]

contact = {692,695,696,698,699,715,719,722,725,746,747,750,753,754}
P = 0.1
sig = []
with open("meme_sites.csv","w",newline="") as f:
    w = csv.writer(f)
    w.writerow(["site","p_value","episodic_positive","is_contact"])
    for i, r in enumerate(rows, start=1):
        p = r[ip]
        e = p <= P
        if e: sig.append(i)
        w.writerow([i, round(p,4), e, i in contact])

print("episodic positive (p<=0.1):", sig)
print("contact among them:", [i for i in sorted(contact) if i in set(sig)] or "none")
