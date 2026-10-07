"""Generate the ring (pyranose / furanose) molecule-library entries for the sugars whose open chain is in js/chem.js.

    python tools/sugar_rings.py            # prints library lines (name|smiles|key|formula|tags|alt|cas); key is left as @KEY@

alpha / beta follow the Haworth convention: ring numbered clockwise from above, a D sugar's anomeric O is DOWN in alpha and UP in beta
(L sugars: the mirror image). Checked against alpha/beta-D-glucopyranose in --check.
"""
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
import numpy as np
RDLogger.DisableLog("rdApp.*")

OPEN = {  # library name -> open-chain SMILES (read from js/chem.js)
}
import re, os
src = open(os.path.join(os.path.dirname(__file__), "..", "js", "chem.js"), encoding="utf-8").read()
for n in "erythrose threose arabinose lyxose ribose xylose allose altrose galactose glucose gulose idose mannose talose fructose fucose rhamnose".split():
    OPEN[n] = re.search(r"^" + n + r"\|([^|]+)\|", src, re.M).group(1)

def chain(m):
    """the carbon backbone from the carbonyl carbon, as atom indices"""
    pat = Chem.MolFromSmarts("[CX3](=O)")
    cands = [a[0] for a in m.GetSubstructMatches(pat)]
    c = [i for i in cands if sum(1 for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol() == "C") <= 2]
    # aldose: carbonyl carbon is a chain end; ketose: C2.  Walk the longest carbon path through it.
    best = []
    def dfs(path):
        nonlocal best
        ext = [n.GetIdx() for n in m.GetAtomWithIdx(path[-1]).GetNeighbors() if n.GetSymbol() == "C" and n.GetIdx() not in path]
        if not ext and len(path) > len(best): best = path[:]
        for e in ext: dfs(path + [e])
    for s in [i for i in range(m.GetNumAtoms()) if m.GetAtomWithIdx(i).GetSymbol() == "C" and m.GetAtomWithIdx(i).GetDegree() <= 2 or False]:
        dfs([s])
    return best

def is_D(m, path):
    """configuration of the highest-numbered stereocentre, by comparing to D-glyceraldehyde"""
    ref = [i for i in path if m.GetAtomWithIdx(i).HasProp("_CIPCode") or True]
    # last carbon bearing OH that is a stereocentre
    sc = [i for i, _ in Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)]
    r = max(sc, key=lambda i: path.index(i))
    q = Chem.MolFromSmarts("[#6][C@H]([OH])[CH2,CH3]")
    q2 = Chem.MolFromSmarts("[#6][C@@H]([OH])[CH2,CH3]")
    hits = [h[1] for h in m.GetSubstructMatches(q, useChirality=True)]
    return r in hits

def ring_forms(open_smiles, size):
    m = Chem.MolFromSmiles(open_smiles)
    path = chain(m)
    if path[0] != [i for i in path][0]: pass
    # make sure the path starts at the carbonyl end
    carbonyl = [a[0] for a in m.GetSubstructMatches(Chem.MolFromSmarts("[CX3]=O"))][0]
    if path[0] != carbonyl and path[-1] == carbonyl: path = path[::-1]
    if carbonyl not in path: raise ValueError("carbonyl not on the chain")
    if path[0] != carbonyl:  # ketose: the carbonyl is C2; start the numbering at the CH2OH end
        pass
    D = is_D(m, path)
    # anomeric carbon index in path; ring closes onto the carbon `size-1` bonds away, through O
    ai = path.index(carbonyl)
    ci = ai + (size - 2)                      # pyranose size=6 -> 4 bonds along; furanose size=5 -> 3
    if ci >= len(path): return None
    anom, closer = path[ai - 0], path[ci]
    ring_o = [n.GetIdx() for n in m.GetAtomWithIdx(closer).GetNeighbors() if n.GetSymbol() == "O"]
    if not ring_o: return None
    ring_o = ring_o[0]
    carb_o = [n.GetIdx() for n in m.GetAtomWithIdx(anom).GetNeighbors() if n.GetSymbol() == "O" and m.GetBondBetweenAtoms(anom, n.GetIdx()).GetBondTypeAsDouble() == 2][0]
    rw = Chem.RWMol(m)
    rw.GetBondBetweenAtoms(anom, carb_o).SetBondType(Chem.BondType.SINGLE)
    rw.AddBond(anom, ring_o, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(carb_o).SetNoImplicit(False); rw.GetAtomWithIdx(ring_o).SetNoImplicit(False)
    Chem.SanitizeMol(rw)
    ring_atoms = path[ai:ci + 1] + [ring_o]
    out = {}
    for tag, name in ((Chem.ChiralType.CHI_TETRAHEDRAL_CW, "cw"), (Chem.ChiralType.CHI_TETRAHEDRAL_CCW, "ccw")):
        x = Chem.Mol(rw); x.GetAtomWithIdx(anom).SetChiralTag(tag); Chem.AssignStereochemistry(x, cleanIt=True, force=True)
        xh = Chem.AddHs(x); AllChem.EmbedMolecule(xh, randomSeed=7); AllChem.MMFFOptimizeMolecule(xh)
        pos = xh.GetConformer().GetPositions()
        P = pos[ring_atoms]; cen = P.mean(0)
        n = sum(np.cross(P[i] - cen, P[(i + 1) % len(P)] - cen) for i in range(len(P)))
        side = {}
        for k, a in enumerate(ring_atoms[:-1]):
            sub = [nb.GetIdx() for nb in xh.GetAtomWithIdx(a).GetNeighbors() if nb.GetIdx() not in ring_atoms and nb.GetSymbol() != "H"]
            if a == anom: sub = [carb_o]
            if sub: side[k + 1] = "down" if np.dot(pos[sub[0]] - pos[a], n) > 0 else "up"
        # normal points away from a viewer looking at a counter-clockwise ring => 'down' for the clockwise Haworth orientation
        out[name] = (Chem.MolToSmiles(x), side[1])
    sx = Chem.MolFromSmiles(open_smiles)
    return D, out, rw, anom

def entries():
    rows = []
    for nm, osm in OPEN.items():
        for size, suf in ((6, "pyranose"), (5, "furanose")):
            r = ring_forms(osm, size)
            if r is None: continue
            D, out, rw, anom = r
            stem = nm[:-4] + "ose" if False else nm
            base = {"fucose": "fuco", "rhamnose": "rhamno"}.get(nm, nm[:-2] if nm.endswith("ose") else nm)
            pre = "D-" if D else "L-"
            alpha_side = "down" if D else "up"
            smi = {k: v for k, v in out.items()}
            a = next(v[0] for v in smi.values() if v[1] == alpha_side)
            b = next(v[0] for v in smi.values() if v[1] != alpha_side)
            un = Chem.Mol(rw); Chem.AssignStereochemistry(un, cleanIt=True, force=True)
            # unspecified anomer: same molecule with the anomeric tag removed
            u = Chem.Mol(rw); u.GetAtomWithIdx(anom).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED); Chem.AssignStereochemistry(u, cleanIt=True, force=True)
            f = rdMolDescriptors.CalcMolFormula(rw)
            full = base + suf
            rows.append((nm, size, pre + full, Chem.MolToSmiles(u), [full], f))
            rows.append((nm, size, "alpha-" + pre + full, a, ["α-" + pre + full, "a-" + pre + full], f))
            rows.append((nm, size, "beta-" + pre + full, b, ["β-" + pre + full, "b-" + pre + full], f))
    return rows

if __name__ == "__main__":
    rows = entries()
    if "--check" in sys.argv:
        g = {r[2]: Chem.MolToSmiles(Chem.MolFromSmiles(r[3])) for r in rows}
        ref = {"alpha-D-glucopyranose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
               "beta-D-glucopyranose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O"}
        for k, v in ref.items(): print(k, g[k] == Chem.MolToSmiles(Chem.MolFromSmiles(v)), g[k])
        print(len(rows), "entries")
    else:
        for nm, size, name, smi, alt, f in rows:
            print("|".join([name, smi, "@KEY@", f, "alcohols,ethers,rings,chiral", ";".join(alt) , ""]))
