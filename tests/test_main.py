import pytest
from fastapi.testclient import TestClient
import os
from src.main import app  
client = TestClient(app)

# Testing Post method
def test_create_item():
    response = client.post(
        '/add',
        json={'id': 5, 'smiles': 'c1cc(C)ccc1'} 
    )
    assert response.status_code == 201  # HTTP 201 Created
    data = response.json()  
    assert data['id'] == 5
    assert data['smiles'] == 'c1cc(C)ccc1'

# Usage of xfail
@pytest.mark.xfail(reason='Already in database')
def test_already_in_db():
    response = client.post(
        '/add',
        json={'id': 1, 'smiles': 'c1cc(C)ccc1'} 
    )
    assert response.status_code == 400  
    data = response.json()  
    assert data['id'] == 5
    assert data['smiles'] == 'c1cc(C)ccc1'

def test_invalid_smiles_add():
    response = client.post(
        '/add',
        json={'id': 7, 'smiles': 'invalid'} 
    )
    assert response.status_code == 400 
    data = response.json()  
    assert data['detail']=='Invalid SMILES string'
# Testing get method using id

def test_valid_id():
    response = client.get(
        '/molecule/1'
    )
    assert response.status_code == 200 
    data = response.json()  
    assert data['id']==1
    assert data['smiles']=='c1cc(C)ccc1'
def test_invalid_id():
    response = client.get(
        '/molecule/sgsgsg'
    )
    assert response.status_code == 422 #Id is not integer in this case 
    data = response.json()  
    assert data['detail']=='Molecule not found'
def test_invalid_id():
    response = client.get(
        '/molecule/9'
    )
    assert response.status_code == 404 # Id is not found
    data = response.json()  
    assert data['detail']=='Molecule not found'

# Tesing Put Method
def test_valid_put():
    response=client.put('/molecules/1',
                       json={'id':1,'smiles':'CCO'})
    assert response.status_code==202
    data=response.json()
    assert data['id']==1
    assert data['smiles']=='CCO'
def test_invalid_smiles_put():
    response=client.put('/molecules/1',
                       json={'id':1,'smiles':'Invalid'})
    assert response.status_code==400
    data=response.json()
    assert data['detail']=='Invalid SMILES string'

def test_mol_not_found_put():
     response=client.put('/molecules/999',
                       json={'id':999,'smiles':'CCO'})
     assert response.status_code==404
     data=response.json()
     assert data['detail']=='Molecule not found'
# Testing Delete Method
def test_valid_id_delete():
    response=client.delete('/molecules/1')
    assert response.status_code==200
# @pytest.mark.xfail(reason='ID does not exist')
def test_invalid_id_detele():
    # It is already deleted
    response=client.delete('/molecules/1') 
    assert response.status_code==404
def test_invalid_id_1_detele():
    # It is already deleted
    response=client.delete('/molecules/999') 
    assert response.status_code==404
# Testing get method of all molecules
def test_get_all():
    response=client.get('/molecules')
    assert response.status_code==200
# Testing Substructure search
def test_valid_search():
    response=client.post('/substructure_search',
                         json={
            'mols': ['c1cc(C)ccc1', 'CCO', 'CC(=O)O', 'CC(=O)Oc1ccccc1C(=O)O'],
            'mol': 'c1cc(C)ccc1'
        }
    )
    assert response.status_code==200
    data=response.json()
    assert isinstance(data,list)
    assert 'c1cc(C)ccc1' in data

def test_invalid_search():
     response=client.post('/substructure_search',
                         json={
            'mols': ['c1cc(C)ccc1', 'CCO', 'CC(=O)O', 'CC(=O)Oc1ccccc1C(=O)O'],
            'mol': 'invalid smiles'
        }
    )
     assert response.status_code==400
     data=response.json()
     assert data['detail']=='Invalid query SMILES string'

#Testing File upload
@pytest.mark.xfail(reason='ID is already existing')
def test_invalid_ids_file_upload():
    file_path=os.path.join(os.path.dirname(__file__),'book2.csv')
    with open(file_path,'rb') as file:
        file_content=file.read()
        response=client.post('/uploadFile', files={'file':('book2.csv', file,'text/csv')})
        assert response.status_code==400
        assert isinstance(response.json(),list)


        
def test_vadid_ids_file_upload():
    file_path=os.path.join(os.path.dirname(__file__),'book3.csv')
    with open(file_path,'rb') as file:
        file_content=file.read()
        response=client.post('/uploadFile', files={'file':('book3.csv', file,'text/csv')})
        assert response.status_code==200
        assert isinstance(response.json(),list)
def test_invalid_smiles_file_upload():
    file_path=os.path.join(os.path.dirname(__file__),'book4.csv')
    with open(file_path,'rb') as file:
        file_content=file.read()
        response=client.post('/uploadFile', files={'file':('book4.csv', file,'text/csv')})
        assert response.status_code==400
