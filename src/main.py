# from fastapi import FastAPI, Request, UploadFile, File, HTTPException, status
# from pydantic import BaseModel
# from rdkit import Chem
# from rdkit.Chem import MolFromSmiles
# from typing import Dict
# import csv
# import io
# from os import getenv


# class Molecul(BaseModel):
#     id: int
#     smiles: str


# app = FastAPI()


# # Initialize molecules_db as a dictionary
# molecules_db: Dict[int, Molecul] = {
#     1: Molecul(id=1, smiles="c1cc(C)ccc1"),
#     2: Molecul(id=2, smiles="CCO"),
#     3: Molecul(id=3, smiles="CC(=O)O"),
#     4: Molecul(id=4, smiles="CC(=O)Oc1ccccc1C(=O)O")
# }


# @app.get("/")
# def get_server():
#     return {"server_id": getenv("SERVER_ID", "1")}


# # Add molecule to the database
# @app.post('/add', status_code=status.HTTP_201_CREATED)
# def add_molecule(molecule: Molecul):
#     if molecule.id in molecules_db:
#         raise HTTPException(status_code=400, detail='Id already exists')
#     if not MolFromSmiles(molecule.smiles):
#         raise HTTPException(
#            status_code=400,
#            detail='Invalid SMILES string')

#     # Make sure that SMILES is valid
#     molecules_db[molecule.id] = molecule

#     return molecule


# # Get Molecule by identifier
# @app.get('/molecule/{molecule_id}')
# def get_molecule(molecule_id: int, summary='Get specific molecule'):
#     '''
#     To Get Molecule By ID
#     '''
#     molecule = molecules_db.get(molecule_id)
#     if molecule:
#         return molecule
#     raise HTTPException(
#         status_code=status.HTTP_404_NOT_FOUND,
#         detail='Molecule not found')


# # Update Molecule by its ID


# @app.put('/molecules/{molecule_id}', status_code=status.HTTP_202_ACCEPTED)
# def update_molecule(molecule_id: int, update_mol: Molecul):
#     '''
#     To Update Molecule By ID
#     '''
#     if molecule_id in molecules_db:
#         if MolFromSmiles(update_mol.smiles):
#             molecules_db[molecule_id] = update_mol
#             return update_mol
#         else:
#             raise HTTPException(
#                 status_code=400,
#                 detail='Invalid SMILES string')
#     raise HTTPException(
#         status_code=status.HTTP_404_NOT_FOUND,
#         detail='Molecule not found'
#         )


# # Delete molecule by ID
# @app.delete('/molecules/{molecule_id}')
# def molecule_deletion(molecule_id: int):
#     '''
#     Delete a molecule by its ID
#     '''
#     if molecule_id in molecules_db:
#         deleted = molecules_db.pop(molecule_id)
#         return deleted
#     raise HTTPException(
#         status_code=status.HTTP_404_NOT_FOUND,
#         detail='Molecule not found')

# # List All molecules


# @app.get('/molecules', status_code=status.HTTP_200_OK)
# def return_molecules():
#     return list(molecules_db.values())

# # Substructure Search


# @app.post('/substructure_search', status_code=status.HTTP_200_OK)
# async def substructure_search(request: Request):
#     '''
#     Performing Substructure Search
#     '''
#     body = await request.json()
#     mols = body.get('mols')
#     mol = body.get('mol')
#     if not mols or not mol:
#         raise HTTPException(
#             status_code=status.HTTP_400_BAD_REQUEST,
#             detail='Missing mols or mol in the request body')
#     molecule = Chem.MolFromSmiles(mol)
#     if molecule is None:
#         raise HTTPException(
#             status_code=status.HTTP_400_BAD_REQUEST,
#             detail='Invalid query SMILES string')
#     match = []
#     for smiles in mols:
#         x = Chem.MolFromSmiles(smiles)
#         if x and x.HasSubstructMatch(molecule):
#             match.append(smiles)

#     return match

# # Upload CSV file and add molecules to the database


# @app.post('/uploadFile')
# async def create_upload(file: UploadFile = File(...)):
#     if file.filename.endswith('.csv'):
#         content = await file.read()
#         # Decoding
#         csv_data = io.StringIO(content.decode('utf-8'))
#         csv_reader = csv.reader(csv_data)
#         added_molecules = []
#         for id, smiles in csv_reader:
#             try:
#                 id = int(id)
#                 if Chem.MolFromSmiles(smiles):
#                     if id not in molecules_db:
#                         molecule = Molecul(id=id, smiles=smiles)
#                         molecules_db[id] = molecule
#                         added_molecules.append(molecule)
#                     else:
#                         raise HTTPException(
#                             status_code=status.HTTP_400_BAD_REQUEST,
#                             detail=f'Molecule with ID {id} already exists')
#             except ValueError:
#                 continue
#         if not added_molecules:
#             raise HTTPException(
#                 status_code=status.HTTP_400_BAD_REQUEST,
#                 detail='No valid molecules added')
#         return added_molecules
#     else:
#         raise HTTPException(
#             status_code=status.HTTP_400_BAD_REQUEST,
#             detail='File must be a CSV')


from fastapi import FastAPI, Depends, HTTPException
from sqlalchemy.orm import Session
from src import models, schema, database

from src.upload import router as upload_router
from src.subs import router as substructure_router


app = FastAPI()
app.include_router(upload_router, prefix="/api", tags=["upload"])
app.include_router(substructure_router, prefix="/api", tags=["substructure"])


models.Base.metadata.create_all(bind=database.engine)


@app.post("/molecules/", response_model=schema.Molecule)
def create_molecule(molecule: schema.CreateMolecule, db: Session = Depends(database.get_db)):
    db_molecule = models.Molecule(name=molecule.name, smiles=molecule.smiles)
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)
    return db_molecule

@app.put('molecules/{molecule.id}',response_model=schema.Molecule)
def update_mol(molecule_id:int,molecule:schema.UpdateMolecule,db:Session=Depends(database.get_db)):
    db_mol=db.query(models.Molecule).filter(models.Molecule.id == molecule_id).first()
    if not db_mol:
        raise HTTPException(status_code=404, detail='Molecule not found')
    db_mol.name=molecule.name
    db_mol.smiles=molecule.smiles
    db.commit()
    db.refresh(db_mol)
    return db_mol

@app.get('molecules/{molecule.id}',response_model=schema.Molecule)
def get_molecule(molecule_id:int, db:Session=Depends(database.get_db)):
    db_mol=db.query(models.Molecule).filter(models.Molecule.id==molecule_id).first()
    if not db_mol:
        raise HTTPException(status_code=404, detail='molecule not found')
    return db_mol
