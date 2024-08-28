# from fastapi import APIRouter, Depends, HTTPException
# from sqlalchemy.orm import Session
# import schema
# from database import get_db

# router = APIRouter()


# @router.post("/molecules/", response_model=schema.Molecule)
# def create_molecule(molecule: schema.MoleculeCreate, db: Session
#                     = Depends(get_db)):
#     return crud.create_molecule(db=db, molecule=molecule)


# @router.get("/molecules/{molecule_id}", response_model=schema.Molecule)
# def read_molecule(molecule_id: int, db: Session = Depends(get_db)):
#     db_molecule = crud.get_molecule(db, molecule_id=molecule_id)
#     if db_molecule is None:
#         raise HTTPException(status_code=404, detail="Molecule not found")
#     return db_molecule
