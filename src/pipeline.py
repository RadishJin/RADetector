# Library
from __future__ import annotations

# Classes
from chembl_webresource_client.new_client import QuerySet

# Custom Modules
from io_parser import io_prot_bindmol, parse_prot_query, io_mol_property, parse_mol_property
from calc import calc_vdwvolume

# Custom Errors
from config import NotFoundError, NotValidError, EmbedFailError


# 전역변수
# =========================================
PROTEIN_IDS : list[str] = [
    'CHEMBL203',    # EFGR
    'CHEMBL1824',   # HER2
    'CHEMBL1844',   # CSF1R
    'CHEMBL301',    # CDK2
    'CHEMBL205'     # CA2
]

OMIT_PROT : list[str] = []
OMIT_MOL : dict[str, list] = {}
RES_DICT : dict = {}
# =========================================


# Pipeline
# =========================================
# for prot_id in PROTEIN_IDS:
#     _omit_mol : list = []

#     # Input Binding Molecules
#     # =====================================
#     try:
#         mol_raw : QuerySet = io_prot_bindmol(prot_id)

#     except NotFoundError as e:
#         OMIT_PROT.append(prot_id)
#         print(f'Proein : {prot_id} : {e}')
#         continue
#     # =====================================


#     # Parse & Create Molecule Property Dictionary
#     # =====================================
#     for raw in mol_raw:
#         mol : tuple = parse_prot_query(raw)


#         # Set variables
#         # =================================
#         mol_id : str = mol[0]
#         prop_dict : dict[str, float] = mol[1]
#         # =================================


#         # Input Molecule Properties
#         # =================================
#         try:
#             prop_raw : QuerySet = io_mol_property(mol_id)

#         except NotFoundError as e:
#             _omit_mol.append(mol_id)
#             print(f'Protein : {prot_id} - Molecule : {mol_id} : {e}')
#             continue
#         # =================================


#         # Parse Fundamental Molecule Properties    
#         # =================================
#         try:
#             prop : dict = parse_mol_property(prop_raw)

#         except (NotFoundError, NotValidError) as e:
#             _omit_mol.append(mol_id)
#             print(f'Protein : {prot_id} - Molecule : {mol_id} : {e}')
#             continue

#         prop_dict.update(prop)
#         # =================================


#         # Calculate Additional Properties
#         # =================================
#         canonical_smiles : str = prop_dict['canonical_smiles']

#         # Van Der Waals Volume
#         try:
#             vdw : float = calc_vdwvolume(canonical_smiles)

#         except (NotValidError, EmbedFailError) as e:
#             _omit_mol.append(mol_id)
#             print(f'Protein : {prot_id} - Molecule : {mol_id} : {e}')
#             continue

#         # Acidity, Basicity
#         # 미완
#         # =================================

#     OMIT_MOL[prot_id] = _omit_mol
#     # =====================================

# =========================================



# Test
# =========================================
test_id = 'CHEMBL203'

raw = io_prot_bindmol(test_id)          # Protein ID 해당하는 모든 Binding Molecules Input - QuerySet >> 여기서 NotFound 발생 가능성
                                                                                                 
molecule = raw[129]                     # 아무 Molecule이나 뽑아서

mol = parse_prot_query(molecule)        # ID, dict(Properties) - Tuple

mol_raw = io_mol_property(mol[0])       # Molecule ID 해당하는 모든 Chemical Properties Input - Dict  >> 여기서 NotFound 발생 가능성

mol_prop = parse_mol_property(mol_raw)  # 필요한 Property만 모아서                                     >> 여기서 NotValid 발생 가능성

mol[1].update(mol_prop)                 # Molecule Tuple의 Property Dict에 추가

mol[1]['vdw'] = calc_vdwvolume(mol[1]['canonical_smiles']) # Van der Waals Volume 계산 후 더해주기     >> 여기서 NotValid, EmbedFail 발생 가능성

print(mol)
# =========================================
