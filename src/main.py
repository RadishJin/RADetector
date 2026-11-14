#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries

# Third-Party Libraries
from chembl_webresource_client.new_client import QuerySet

# Custom Modules & Classes
from io_parser import io_prot_bindmol, parse_mol_query
from pipeline import get_known_props
from build import how_similar
from exceptions import PipelineError
from config import Setting, Config
from domain import Molecule, RADetect

#endregion
#region MAIN

setting = Setting(target_proteins= ['CHEMBL203'])
config = Config(setting.target_proteins)

# for prot_id in config.prot_ids:

#     try:
#         mol_raw : QuerySet = io_prot_bindmol(prot_id)
#     except PipelineError as e:
#         config.omits[prot_id] = 'All'
#         print(f'Protein : {prot_id} : {e}')
#         continue 

#     config.omits[prot_id] = []
#     for raw in mol_raw:
#         data : tuple = parse_mol_query(raw)

#         try:
#             mol = Molecule.extended(data[0], data[1], data[2], data[3], setting.maxseed, setting.maxiter)

#         except PipelineError as e:
#                     config.omits[prot_id].append(data[0])
#                     print(f'Protein : {prot_id} - Molecule : {data[0]} : {e}')
#                     continue
        
#         try:
#             mol.properties.update(get_known_props(mol.id))

#         except PipelineError as e:
#             config.omits[prot_id].append(mol.id)
#             print(f'Protein : {prot_id} - Molecule : {mol.id} : {e}')
#             continue

#         print(f'Successful with {prot_id} - {mol.id}')
#         config.res_mols.append(mol)
    
#         # 아래로는 미정

# print(len(mol_raw))
# print(len(config.omits))
# print(len(config.res_mols))

#endregion

# test
prot_id = config.prot_ids[0]

mol_raw : QuerySet = io_prot_bindmol(prot_id)

for i in range(10, 224):
    data : tuple = parse_mol_query(mol_raw[i])
    mol = Molecule.extended(data[0], data[1], data[2], data[3], setting.maxseed, setting.maxiter)
    mol.properties.update(get_known_props(mol.id))
    config.res_mols.append(mol)

rad_a = RADetect.extended(prot_id, config.res_mols)


prot_a = rad_a.array

print(prot_a)
print(how_similar(prot_a, prot_a, gamma = 1.0))






# setting에서 프로틴 아이디를 받게하기, config는 몰라야함
# batch API
# MultiProcessing