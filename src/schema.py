from pydantic import BaseModel


class CreateMolecule(BaseModel):
    name: str
    smiles: str
    # formula: str


class UpdateMolecule(BaseModel):
    name: str
    smiles: str
    # formula: str


class Molecule(BaseModel):
    id: int
    name: str
    smiles: str
    # formula: str

    class Config:
        orm_mode = True
