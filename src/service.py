from typing import TYPE_CHECKING, List

from .database import database
from .models import models
from .schema import schema

if TYPE_CHECKING:
    from sqlalchemy.orm import Session


def _add_tables():
    return database.Base.metadata.create_all(bind=database.engine)


def get_db():
    db = database.SessionLocal()
    try:
        yield db
    finally:
        db.close()


async def create_molecule(
    molecule: schema.CreateMolecule, db: "Session"
) -> schema.Molecule:
    mol = models.Molecule(**molecule.dict())
    db.add(mol)
    db.commit()
    db.refresh(mol)
    return schema.Molecule.from_orm(mol)


async def get_all_molecules(db: "Session") -> List[schema.Molecule]:
    molecules = db.query(models.Molecule).all()
    return list(map(schema.Molecule.from_orm, molecules))


async def get_molecule(molecule_id: int, db: "Session") -> schema.Molecule:
    molecule = db.query(models.Molecule).filter(
        models.Molecule.id == molecule_id).first()
    if molecule:
        return schema.Molecule.from_orm(molecule)
    return None


async def delete_molecule(molecule: models.Molecule, db: "Session"):
    db.delete(molecule)
    db.commit()


async def update_molecule(
    molecule_data: schema.UpdateMolecule,
    molecule: models.Molecule, db: "Session"
) -> schema.Molecule:
    molecule.name = molecule_data.name
    molecule.smiles = molecule_data.smiles
    molecule.formula = molecule_data.formula

    db.commit()
    db.refresh(molecule)

    return schema.Molecule.from_orm(molecule)
