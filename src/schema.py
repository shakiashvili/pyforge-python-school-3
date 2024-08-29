from pydantic import BaseModel


class CreateMolecule(BaseModel):
    name: str
    smiles: str


class UpdateMolecule(BaseModel):
    name: str
    smiles: str


class Molecule(BaseModel):
    id: int
    smiles: str

    class Config:
        orm_mode = True