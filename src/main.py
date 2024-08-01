from fastapi import FastAPI, UploadFile, File, HTTPException, status
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import Draw
from typing import Dict, List
import csv
import io
from os import getenv

class Molecul(BaseModel):
    id: int
    smiles: str

app = FastAPI()

# Initialize molecules_db as a dictionary
molecules_db: Dict[int, Molecul] = {
    1: Molecul(id=1, smiles="c1cc(C)ccc1"),
    2: Molecul(id=2, smiles="CCO"),
    3: Molecul(id=3, smiles="CC(=O)O"),
    4: Molecul(id=4, smiles="CC(=O)Oc1ccccc1C(=O)O")
}

@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}

# Add molecule to the database
@app.post('/add', status_code=status.HTTP_201_CREATED)
def add_molecule(molecule: Molecul):
    if molecule.id in molecules_db:
        raise HTTPException(status_code=400, detail='Id already exists')
    molecules_db[molecule.id] = molecule
    return molecule

# Get Molecule by identifier
@app.get('/molecule/{molecule_id}')
def get_molecule(molecule_id: int, summary='Get specific molecule'):
    '''
    To Get Molecule By ID
    '''
    molecule = molecules_db.get(molecule_id)
    if molecule:
        return molecule
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail='Molecule not found')

# Update Molecule by its ID
@app.put('/molecules/{molecule_id}', status_code=status.HTTP_202_ACCEPTED)
def update_molecule(molecule_id: int, update_mol: Molecul):
    '''
    To Update Molecule By ID
    '''
    if molecule_id in molecules_db:
        molecules_db[molecule_id] = update_mol
        return update_mol
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail='Molecule not found')

# Delete molecule by ID
@app.delete('/molecules/{molecule_id}')
def molecule_deletion(molecule_id: int):
    '''
    Delete a molecule by its ID
    '''
    if molecule_id in molecules_db:
        deleted = molecules_db.pop(molecule_id)
        return deleted
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail='Molecule not found')

# List All molecules
@app.get('/molecules', status_code=status.HTTP_200_OK)
def return_molecules():
    return list(molecules_db.values())

# Substructure Search
@app.post('/substructure_search', status_code=status.HTTP_200_OK)
def substructure_search(smiles: str):
    '''
    Performing Substructure Search
    '''
    try:
        query_mol = Chem.MolFromSmiles(smiles)
        if query_mol is None:
            raise ValueError('Invalid SMILES string')
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    
    matches = []
    for molecule in molecules_db.values():
        target_mol = Chem.MolFromSmiles(molecule.smiles)
        if target_mol and target_mol.HasSubstructMatch(query_mol):
            matches.append(molecule)
    
    return matches

# Upload CSV file and add molecules to the database
@app.post('/uploadFile')
async def create_upload(file: UploadFile = File(...)):
    if file.filename.endswith('.csv'):
        content = await file.read()
        # Decoding
        csv_data = io.StringIO(content.decode('utf-8'))
        csv_reader = csv.reader(csv_data)
        
        # Skipping the header
        header = next(csv_reader, None)
        
        added_molecules = []
        for id, smiles in csv_reader:
            try:
                id = int(id)
                if Chem.MolFromSmiles(smiles):
                    if id not in molecules_db:
                        molecule = Molecul(id=id, smiles=smiles)
                        molecules_db[id] = molecule
                        added_molecules.append(molecule)
                    else:
                        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f'Molecule with ID {id} already exists')
            except ValueError:
                continue
        
        if not added_molecules:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail='No valid molecules added')
        
        return added_molecules
    else:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail='File must be a CSV')