from fastapi import APIRouter, UploadFile, File, HTTPException, status, Depends
from sqlalchemy.orm import Session
import csv
import io
from rdkit import Chem
from src.database import get_db
from src.schema import CreateMolecule
from src.models import Molecule

router = APIRouter()

@router.post("/uploadFile", status_code=status.HTTP_201_CREATED)
async def create_upload(file: UploadFile = File(...), db: Session = Depends(get_db)):
    if not file.filename.endswith('.csv'):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='File must be a CSV'
        )

    content = await file.read()
    csv_data = io.StringIO(content.decode('utf-8'))
    csv_reader = csv.reader(csv_data)
    
    headers = next(csv_reader, None)  # Skip header row if present
    if headers is None or len(headers) < 2:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='CSV file must contain at least two columns (id and smiles)'
        )
    
    added_molecules = []
    for row in csv_reader:
        if len(row) < 2:
            continue  # Skip rows that don't have enough columns

        try:
            id = int(row[0])
            smiles = row[1]
            if Chem.MolFromSmiles(smiles):
                molecule_data = MoleculeCreate(id=id, smiles=smiles)
                db_molecule = Molecule(**molecule_data.dict())
                db.add(db_molecule)
                added_molecules.append(molecule_data)
        except ValueError:
            continue

    if not added_molecules:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='No valid molecules added'
        )
    
    db.commit()
    return added_molecules