#!/usr/bin/env python3
"""
compute_ddg.py
--------------
Per-residue binding ddG (DDG) scan of the PB1-PB2 interface with EvoEF2.

Works for either side of the interface of PDB 3A1G:
  chain A = PB1 (C-terminal domain)   chain B = PB2 (N-terminal domain)

For each interface contact residue of the chosen chain, mutate it to all 19
other amino acids, rebuild the complex, and compute the
change in PB1-PB2 binding free energy relative to wild type:

    ddG_bind(mut) = E_bind(mutant) - E_bind(wild_type)

Positive ddG_bind  -> mutation DESTABILIZES the interface (harder to tolerate)
Negative / ~0      -> mutation is tolerated at the binding level

Output: CSV with columns  pos, wt, mut, ddG_bind   (n_residues x 19 rows)

Requires:
  - EvoEF2 binary (build once with run_ddg.sh)
  - A repaired structure: 3A1G_clean_Repair.pdb  (produced by run_ddg.sh)

Usage (normally called by run_ddg.sh):
  PB1: python3 compute_ddg.py --evoef EvoEF2/EvoEF2 \
         --pdb struct/3A1G_clean_Repair.pdb --chain A --out ddg_all.csv
  PB2: python3 compute_ddg.py --evoef EvoEF2/EvoEF2 \
         --pdb struct/3A1G_clean_Repair.pdb --chain B --out ddg_all_PB2.csv

Which residues are scanned (--residues):
  default     (default) chain A -> the fixed list of 41 PB1 residues below
              (what the paper used); chain B -> auto.
  auto        every residue of --chain with any heavy atom within --cutoff
              (4.5 A) of the partner chain. NOTE: for chain A this gives 39
              residues, not exactly the 41 of the fixed list (it misses
              741, 745, 749, 752 and adds 729, 735; the original list used a
              slightly different contact definition). Use the fixed list
              for PB1 so results stay comparable.
  "1,5,9"     an explicit comma-separated list of residue numbers.
"""
import argparse, os, re, subprocess, sys, tempfile, csv

# PB1 interface contact residues: 4.5 A to native PB2 (same as --residues auto, chain A)
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

def interface_residues(pdb, chain, cutoff=4.5):
    """Residue numbers of `chain` with any heavy atom within `cutoff` A of any
    heavy atom of the other chain(s). Pure python, no dependencies."""
    mine, other = {}, []
    for line in open(pdb):
        if not line.startswith("ATOM"): continue
        if line[12:16].strip().startswith("H"): continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if line[21] == chain:
            mine.setdefault(int(line[22:26]), []).append(xyz)
        else:
            other.append(xyz)
    c2 = cutoff * cutoff
    hits = []
    for resi, atoms in mine.items():
        if any((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 <= c2
               for a in atoms for b in other):
            hits.append(resi)
    return sorted(hits)

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
    ap.add_argument("--chain", default="A", help="chain to mutate: A = PB1, B = PB2")
    ap.add_argument("--residues", default="default",
                    help="default | auto | comma-separated residue numbers")
    ap.add_argument("--cutoff", type=float, default=4.5,
                    help="contact distance (A) for --residues auto")
    ap.add_argument("--limit", type=int, default=0,
                    help="scan only the first N residues (quick smoke test)")
    ap.add_argument("--out", default="ddg_all.csv")
    a = ap.parse_args()

    evoef = os.path.abspath(a.evoef); pdb = os.path.abspath(a.pdb)
    seq = read_chain_seq(pdb, a.chain)
    wt_bind = binding_energy(evoef, pdb)
    print(f"WT binding dG = {wt_bind:.2f}", file=sys.stderr)

    if a.residues == "default":
        a.residues = "pb1-fixed" if a.chain == "A" else "auto"
    if a.residues == "auto":
        # keep native numbering only (3A1G chain B starts with expression-tag residues <= 0)
        residues = [r for r in interface_residues(pdb, a.chain, a.cutoff) if r >= 1]
    elif a.residues == "pb1-fixed":
        residues = INTERFACE
    else:
        residues = [int(x) for x in a.residues.split(",")]
    if a.limit:
        residues = residues[:a.limit]
    print(f"chain {a.chain}: scanning {len(residues)} residues: {residues}", file=sys.stderr)

    rows = []
    for pos in residues:
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
