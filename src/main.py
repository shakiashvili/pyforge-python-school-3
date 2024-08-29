# from fastapi import FastAPI
# from api import molecules as _molecules
# from api import search as _search
# from api import upload as _upload
# from database import engine as _engine
# from .models import Base as _Base
# _Base.metadata.create_all(bind=_engine)

# app = FastAPI()

# app.include_router(_molecules.router)
# app.include_router(_search.router)
# app.include_router(_upload.router)

from fastapi import FastAPI, Depends, HTTPException
from sqlalchemy.orm import Session
from typing import List
from service import (
    create_molecule,
    get_all_molecules,
    get_molecule,
    delete_molecule,
    update_molecule
)
import schema as _schemas
import database as _database

app = FastAPI()

# Dependency


def get_db():
    db = _database.SessionLocal()
    try:
        yield db
    finally:
        db.close()


@app.post("/molecules/", response_model=_schemas.Molecule)
async def create_molecule_endpoint(
    molecule: _schemas.CreateMolecule, db: Session = Depends(get_db)
):
    return await create_molecule(molecule, db)


@app.get("/molecules/", response_model=List[_schemas.Molecule])
async def read_molecules(db: Session = Depends(get_db)):
    return await get_all_molecules(db)


@app.get("/molecules/{molecule_id}", response_model=_schemas.Molecule)
async def read_molecule(molecule_id: int, db: Session = Depends(get_db)):
    molecule = await get_molecule(molecule_id, db)
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule


@app.delete("/molecules/{molecule_id}", response_model=_schemas.Molecule)
async def delete_molecule_endpoint(molecule_id: int, db: Session
                                   = Depends(get_db)):
    molecule = await get_molecule(molecule_id, db)
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    await delete_molecule(molecule, db)
    return molecule


@app.put("/molecules/{molecule_id}", response_model=_schemas.Molecule)
async def update_molecule_endpoint(
    molecule_id: int, molecule_data: _schemas.UpdateMolecule,
    db: Session = Depends(get_db)
):
    molecule = await get_molecule(molecule_id, db)
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return await update_molecule(molecule_data, molecule, db)
