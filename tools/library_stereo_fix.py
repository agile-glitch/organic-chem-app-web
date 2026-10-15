"""One-off correction of library compounds whose stereo was wrong or missing (found by comparing the library with
ChEMBL's approved drugs, tools/library_add.py). Each new SMILES is ChEMBL's (or PubChem's for LSD), checked here by its
full InChIKey. The library key (stereo-blind) and the formula do not change, so only the SMILES field is rewritten.
  python tools/library_stereo_fix.py [--write]
"""
import os, sys
from rdkit import Chem
from rdkit.Chem import inchi, rdMolDescriptors
HERE = os.path.dirname(os.path.abspath(__file__)); APP = os.path.dirname(HERE)
JS = os.path.join(APP, 'js', 'chem.js')
# library name: (source, its name there, the correct full InChIKey)
FIX = {
    'lsd': ('pubchem', 'CCN(CC)C(=O)[C@H]1CN([C@@H]2CC3=CNC4=CC=CC(=C34)C2=C1)C', 'VAYOSLLFUXYJDT-RDTXWAMCSA-N'),   # (6aR,9R); was its mirror image
    'adrenaline': ('chembl', 'epinephrine', 'UCTWMZQNUQWSLP-VIFPVBQESA-N'),             # (R); stereo was missing
    'nicotine': ('chembl', 'nicotine', 'SNICXCGAKADSCV-JTQLQIEISA-N'),                  # (S); missing
    'tartaric acid': ('chembl', 'tartaric acid', 'FEWJPZIEWOKRBE-LWMBPPNESA-N'),        # L-(+), (2R,3R); missing
    'cocaine': ('chembl', 'cocaine', 'ZPUCINDJVBIVPJ-LJISPDSOSA-N'),                    # (1R,2R,3S,5S); was a diastereomer
    'levetiracetam': ('chembl', 'levetiracetam', 'HPHUVLMMVZITSG-LURJTMIESA-N'),        # (S); was (R)
    'methamphetamine': ('chembl', 'methamphetamine', 'MYWUZJCMWCOHBA-VIFPVBQESA-N'),    # (S); was (R)
    'ephedrine': ('chembl', 'ephedrine', 'KWGRBVOPPLSCSI-WPRPVWTQSA-N'),                # (1R,2S); was (1R,2R)
    'pseudoephedrine': ('chembl', 'pseudoephedrine', 'KWGRBVOPPLSCSI-WCBMZHEXSA-N'),    # (1S,2S); was (1R,2S) = ephedrine
    'atropine': ('chembl', 'atropine', 'RKUNBYITZUJHSG-SPUOUPEWSA-N'),                  # tropane endo ester; missing
    'scopolamine': ('chembl', 'scopolamine', 'STECJAGHUSJQJN-FWXGHANASA-N'),            # missing
    'delta-9-tetrahydrocannabinol': ('chembl', 'dronabinol', 'CYQFCXCEBYINGO-IAGOWNOFSA-N'),   # (6aR,10aR); was another stereoisomer
    'tetracycline': ('chembl', 'tetracycline', 'OFVLGDICTFRJMM-WESIUVDSSA-N'),          # missing
    'spironolactone': ('chembl', 'spironolactone', 'LXMSZDCAJNLERA-ZHYRCANASA-N'),      # C8 was inverted
    'tadalafil': ('chembl', 'tadalafil', 'WOXKDUGGOYFFRN-IIBYNOLFSA-N'),                # (6R,12aR); was (6S,12aR)
    'esomeprazole': ('chembl', 'esomeprazole', 'SUBDBMMJDZJVOS-DEOSSOPVSA-N'),          # (S); was (R)
    'montelukast': ('chembl', 'montelukast', 'UCHDWCPVSPXUMX-TZIWLTJVSA-N'),            # (R); was (S)
    'podophyllotoxin': ('chembl', 'podofilox', 'YJGVMLPVUAXIQN-XVVDYKMHSA-N'),          # was another stereoisomer
    'cholecalciferol': ('chembl', 'cholecalciferol', 'QYSXJUFSXHHAJI-YRZJJWOYSA-N'),    # the CD ring was wrong
}


def main():
    import pandas as pd, datamol
    d = pd.read_parquet(os.path.join(os.path.dirname(datamol.__file__), 'data', 'chembl_approved_drugs.parquet'))
    chembl = {r.pref_name.lower(): r.smiles for _, r in d.iterrows()}
    js = open(JS, encoding='utf-8').read()
    start = js.index('const MOL_LIBRARY_RAW = `') + 25; end = js.index('`;', start)
    lines = js[start:end].split('\n'); done = 0
    for i, l in enumerate(lines):
        f = l.split('|')
        if f[0] not in FIX: continue
        src, ref, want = FIX[f[0]]
        smi = chembl[ref] if src == 'chembl' else ref
        m = Chem.MolFromSmiles(smi)
        if f[0] == 'esomeprazole':          # the library writes sulfoxides S=O, ChEMBL [S+][O-]: same InChI
            m = Chem.MolFromSmiles(smi.replace('[S@@+]([O-])', '[S@@](=O)'))
        assert inchi.MolToInchiKey(m) == want, (f[0], inchi.MolToInchiKey(m))
        assert rdMolDescriptors.CalcMolFormula(m) == f[3], (f[0], f[3])
        Chem.Kekulize(m, clearAromaticFlags=True)
        new = Chem.MolToSmiles(m, kekuleSmiles=True)
        print(f'{f[0]:30s} {f[1]}\n{"":30s} {new}')
        f[1] = new; lines[i] = '|'.join(f); done += 1
    assert done == len(FIX), done
    if '--write' in sys.argv:
        open(JS, 'w', encoding='utf-8').write(js[:start] + '\n'.join(lines) + js[end:]); print(f'{done} entries corrected in js/chem.js')


if __name__ == '__main__':
    main()
