# from fastapi import APIRouter, HTTPException, Request, status
# from fastapi import Depends
# from rdkit import Chem
# from sqlalchemy.orm import Session
# from database import get_db

# router = APIRouter()


# @router.post("/substructure_search", status_code=status.HTTP_200_OK)
# async def substructure_search(request: Request, db: Session = Depends(get_db)):
#     body = await request.json()
#     mols = body.get('mols')
#     mol = body.get('mol')
#     if not mols or not mol:
#         raise HTTPException(
#             status_code=status.HTTP_400_BAD_REQUEST,
#             detail='Missing mols or mol in the request body')
#     query_mol = Chem.MolFromSmiles(mol)
#     if query_mol is None:
#         raise HTTPException(
#             status_code=status.HTTP_400_BAD_REQUEST,
#             detail='Invalid query SMILES string')
#     matches = []
#     for smiles in mols:
#         molecule = Chem.MolFromSmiles(smiles)
#         if molecule and molecule.HasSubstructMatch(query_mol):
#             matches.append(smiles)

#     return matches
