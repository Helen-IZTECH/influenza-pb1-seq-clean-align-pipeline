#!/usr/bin/env python3
"""
compute_ddg.py
--------------
Per-residue binding ddG (DDG) scan of the PB1-PB2 interface with EvoEF2.

For each of the 41 PB1 interface contact residues (chain A of PDB 3A1G),
mutate it to all 19 other amino acids, rebuild the complex, and compute the
change in PB1-PB2 binding free energy relative to wild type:

    ddG_bind(mut) = E_bind(mutant) - E_bind(wild_type)

Positive ddG_bind  -> mutation DESTABILIZES the interface (harder to tolerate)
Negative / ~0      -> mutation is tolerated at the binding level

Output: ddg_all.csv with columns  pos, wt, mut, ddG_bind   (41 x 19 = 779 rows)

Requires:
  - EvoEF2 binary (build once with run_ddg.sh)
  - A repaired structure: 3A1G_clean_Repair.pdb  (produced by run_ddg.sh)

Usage:
  python3 compute_ddg.py --evoef EvoEF2/EvoEF2 \
      --pdb struct/3A1G_clean_Repair.pdb --chain A --out ddg_all.csv
"""
import argparse, os, re, subprocess, sys, tempfile, csv

# PB1 interface contact residues: 4.5 A to native PB2, from find_interface_contacts.R
INTERFACE = [685,686,688,689,690,691,692,693,695,696,697,698,699,708,709,710,
             711,712,713,714,715,717,718,719,721,722,725,739,741,742,743,745,
             746,747,749,750,752,753,754,756,757]

AA3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
          'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
          'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
AA20 = list("ACDEFGHIKLMNPQRSTVWY")

BIND_RE = re.compile(r"Total\s+binding.*?=\s*(-?\d+\.\d+)", re.IGNORECASE)

def read_chain_seq(pdb, chain):
    """Map residue number -> 1-letter code for CA atoms of `chain`."""
    seq = {}
    for line in open(pdb):
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
            resn = line[17:20].strip(); resi = int(line[22:26])
            if resn in AA3to1:
                seq[resi] = AA3to1[resn]
    return seq

def binding_energy(evoef, pdb):
    """Run EvoEF2 ComputeBinding on a 2-chain complex, return total binding dG."""
    out = subprocess.run([evoef, "--command=ComputeBinding", "--pdb="+pdb],
                         capture_output=True, text=True).stdout
    m = BIND_RE.search(out)
    if not m:  # fall back: EvoEF2 prints "Total  = ..." for the binding block
        m = re.search(r"Total\s*=\s*(-?\d+\.\d+)", out)
    return float(m.group(1)) if m else None

def build_mutant(evoef, pdb, mut_tag, workdir):
    """BuildMutant for a single mutation, return path to the rebuilt PDB."""
    indiv = os.path.join(workdir, "individual_list.txt")
    open(indiv, "w").write(mut_tag + ";\n")
    subprocess.run([evoef, "--command=BuildMutant", "--pdb="+pdb,
                    "--mutant_file="+indiv], cwd=workdir,
                   capture_output=True, text=True)
    base = os.path.splitext(os.path.basename(pdb))[0]
    cand = os.path.join(workdir, base + "_Model_0001.pdb")
    return cand if os.path.exists(cand) else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evoef", required=True)
    ap.add_argument("--pdb", required=True, help="repaired complex PDB")
    ap.add_argument("--chain", default="A", help="PB1 chain id")
    ap.add_argument("--out", default="ddg_all.csv")
    a = ap.parse_args()

    evoef = os.path.abspath(a.evoef); pdb = os.path.abspath(a.pdb)
    seq = read_chain_seq(pdb, a.chain)
    wt_bind = binding_energy(evoef, pdb)
    print(f"WT binding dG = {wt_bind:.2f}", file=sys.stderr)

    rows = []
    for pos in INTERFACE:
        wt = seq.get(pos)
        if wt is None:
            print(f"  skip {pos}: not resolved in chain {a.chain}", file=sys.stderr); continue
        for mut in AA20:
            if mut == wt: continue
            tag = f"{wt}{a.chain}{pos}{mut}"   # EvoEF2 mutant syntax: {WT}{Chain}{Pos}{Mut}
            with tempfile.TemporaryDirectory() as wd:
                mpdb = build_mutant(evoef, pdb, tag, wd)
                ddg = ""
                if mpdb:
                    mb = binding_energy(evoef, mpdb)
                    if mb is not None:
                        ddg = round(mb - wt_bind, 2)
                rows.append((pos, wt, mut, ddg))
        print(f"  done {pos} ({wt})", file=sys.stderr)

    with open(a.out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["pos","wt","mut","ddG_bind"]); w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {a.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
