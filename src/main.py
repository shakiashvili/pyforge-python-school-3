from rdkit import Chem
from rdkit.Chem import Draw
from fastapi import FastAPI,UploadFile,status,HTTPException,File
from model import Molecul
from typing import List
import csv
import io
# def substructure_search(mols, mol):
#     molecule = Chem.MolFromSmiles(mol)
#     match = []
#     for smiles in mols:
#         x = Chem.MolFromSmiles(smiles)
#         if x and x.HasSubstructMatch(molecule):
#             match.append(smiles)
#     return match

# mols = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
# mol = "c1ccccc1"
# result=substructure_search(mols, mol)
# print(result)
# molecules=[Chem.MolFromSmiles(smiles) for smiles in result]
# img=Draw.MolsToGridImage(molecules)
# img.show()

app=FastAPI()


molecules_db:List [Molecul] = [
    Molecul(id=1, smiles="c1cc(C)ccc1"),
    Molecul(id=2, smiles="CCO"),
    Molecul(id=3, smiles="CC(=O)O"),
    Molecul(id=4, smiles="CC(=O)Oc1ccccc1C(=O)O")
]

#Add molecule to the database
@app.post('/add', status_code=status.HTTP_201_CREATED)
def add_molecule(molecule: Molecul):
    if any(mol.id==molecule.id for mol in molecules_db):
        raise HTTPException(status_code=400, detail='Id already exists')
    molecules_db.append(molecule.dict())
    return molecule
#Get Molecule by identifier
@app.get('/molecule/{molecule_id}')
def get_molecule(molecule_id:int,summary='Get specific molecule'):
    '''
    To Get Molecule By ID
    '''
    for molecule in molecules_db:
        if molecule.id==molecule_id:
            return molecule
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,detail='Molecule is not found')

#Update Molecule by it's ID
@app.put('/molecules/{molecule_id}')
def update_molecule(molecule_id: int,update_mol:Molecul,status_code=202):
    '''
    To Update Molecule By ID
    '''
    for index,mol in enumerate(molecules_db):
        if mol.id==molecule_id:
            molecules_db[index]=update_mol
            return update_mol
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,detail='Molecule is not found')

#Delete molecule by ID
@app.delete('/molecules/{molecule_id}')
def molecule_deletion(molecule_id:int):
    for index,mol in enumerate(molecules_db):
        if mol.id==molecule_id:
            deleted=molecules_db.pop(index)
            return deleted
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,detail='Molecule not found')
#List All molecules
@app.get('/molecules',status_code=201)
def return_molecules():
    return molecules_db
#Substructure Search for
@app.post('/substucture_search',status_code=status.HTTP_200_OK)
def substructure_search(smiles:str):
    '''
    Perfomring Substructure Search
    '''
    try:
       mol=Chem.MolFromSmiles(smiles)
       if mol is None:
         raise ValueError('Invalid Smile String')
    except Exception as e:
        raise HTTPException('Not Valid Smile string',e)
    match=[]
    for molecule in molecules_db:
        target=Chem.MolFromSmiles(molecule.smiles)
        if target and target.HasSubstructMatch(mol):
            match.append(molecule)
    return match

@app.post('/uploadFile')
async def create_upload(file:UploadFile=File(...)):
    if file.filename.endswith('csv'):
        content=await file.read()
        # Decoding 
        csv_data=io.StringIO(content.decode('Utf-8'))
        csv_reader=csv.DictReader(csv_data)
        # Initializing empty list 
        smil=[]
        for row in csv_reader:
            if 'id' in row and 'smiles' in row:
                try:
                   id=row['id']
                   smiles=row['smiles']
                   if Chem.MolFromSmiles(smiles):
                        molecule=Molecul(id=id,smiles=smiles)
                        molecules_db.append(molecule)
                        smil.append(molecule)
                except ValueError:
                    continue
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST,detail='File must be in specific Format')
    return smil