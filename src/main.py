from rdkit import Chem
from rdkit.Chem import Draw
def substructure_search(mols, mol):
    molecule = Chem.MolFromSmiles(mol)
    match = []
    for smiles in mols:
        x = Chem.MolFromSmiles(smiles)
        if x and x.HasSubstructMatch(molecule):
            match.append(smiles)
    return match

mols = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
mol = "c1ccccc1"
result=substructure_search(mols, mol)
print(result)
molecules=[Chem.MolFromSmiles(smiles) for smiles in result]
img=Draw.MolsToGridImage(molecules)
img.show()