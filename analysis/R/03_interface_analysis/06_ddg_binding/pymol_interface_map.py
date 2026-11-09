#!/usr/bin/env python3
"""
pymol_interface_map.py
----------------------
PyMOL structural map of the PB1-PB2 interface (PDB 3A1G), residues coloured by
(1) selection class from stage 05 and (2) mean binding ddG from stage 06.

Inputs : struct/3A1G_clean.pdb  (made by run_ddg.sh)
         ddg_vs_selection.csv   (made by merge_ddg_selection.R)
Outputs: FIG_pymol_selection.png, FIG_pymol_ddG.png, interface_map.pse

Run (needs PyMOL; `pip install pymol-open-source-whl` or conda `pymol-open-source`):
    python3 pymol_interface_map.py
    python3 pymol_interface_map.py --pdb struct/3A1G_clean.pdb --csv ddg_vs_selection.csv
"""
import argparse, csv
import pymol
from pymol import cmd

ap = argparse.ArgumentParser()
ap.add_argument("--pdb", default="struct/3A1G_clean.pdb")
ap.add_argument("--csv", default="ddg_vs_selection.csv")
ap.add_argument("--chain", default="A", help="chain the CSV refers to (A = PB1, B = PB2)")
ap.add_argument("--highlight", default="715", help="comma list of residues to label")
a = ap.parse_args()

pymol.finish_launching(["pymol", "-cq"])          # headless
cmd.load(a.pdb, "cx")
cmd.remove("hydro or name H*")
partner = "B" if a.chain == "A" else "A"
cmd.hide("everything"); cmd.bg_color("white")
cmd.show("cartoon", "cx")
cmd.color("grey70", f"cx and chain {a.chain}")
cmd.color("lightblue", f"cx and chain {partner}")
cmd.set("cartoon_transparency", 0.5); cmd.set("stick_radius", 0.22)
cmd.set("ray_opaque_background", 1); cmd.set("antialias", 2)

rows = list(csv.DictReader(open(a.csv)))
CLASS_COL = {"invariant": "firebrick", "purifying": "orange", "tolerant": "marine"}
for r in rows:
    sel = f"cx and chain {a.chain} and resi {r['pos']}"
    cmd.show("sticks", sel)
    cmd.color(CLASS_COL.get(r["class"], "grey50"), sel + " and elem C")
for p in a.highlight.split(","):
    cmd.label(f"cx and chain {a.chain} and resi {p} and name CA", '"%s%s" % (resn, resi)')
cmd.set("label_size", 24); cmd.set("label_color", "black")
cmd.orient("cx"); cmd.zoom(f"cx and chain {a.chain} and resi " + "+".join(r["pos"] for r in rows), 8)
cmd.png("FIG_pymol_selection.png", 1600, 1200, dpi=200, ray=1)

# second colouring: mean binding ddG (blue = tolerated, red = destabilising)
cmd.alter(f"cx and chain {a.chain}", "b=0")
for r in rows:
    cmd.alter(f"cx and chain {a.chain} and resi {r['pos']}", "b=%s" % r["mean_ddG"])
cmd.spectrum("b", "blue_white_red", f"cx and chain {a.chain} and resi " + "+".join(r["pos"] for r in rows),
             minimum=-2, maximum=8)
cmd.png("FIG_pymol_ddG.png", 1600, 1200, dpi=200, ray=1)
cmd.save("interface_map.pse")
print("wrote FIG_pymol_selection.png, FIG_pymol_ddG.png, interface_map.pse")
cmd.quit()
