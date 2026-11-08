#region MODULE & LIBRARY

from __future__ import annotations

# Standard Libraries
# ============================
# ============================

# Third-Party Libraries
# ============================
from chembl_webresource_client.new_client import QuerySet
# ============================

# Custom Modules & Classes
# ============================
from io_parser import io_prot_bindmol, parse_prot_query, io_mol_property, parse_mol_property
from calc import calc_vdwvolume

from config import PipelineError
# ============================

#endregion
#region GLOBAL VARIABLE

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

#endregion
#region INTERNAL

# ============================1
# ============================1

#endregion
#region MAIN

# ============================
# # def pipeline():
# for prot_id in PROTEIN_IDS:
#     _omit_mol : list = []

#     # 1단계: Input Binding Molecules
#     # ============================

#     try:
#         mol_raw : QuerySet = io_prot_bindmol(prot_id)

#     except PipelineError as e:
#         OMIT_PROT.append(prot_id)
#         print(f'Protein : {prot_id} : {e}')
#         continue

#     # ============================


#     # 2단계: Parse & Create Molecule Property Dictionary
#     # ============================

#     for raw in mol_raw:

#             # 2-0 단계 : Set Variable
#             mol : tuple = parse_prot_query(raw)
#             mol_id : str = mol[0]
#             prop_dict : dict[str, float] = mol[1]
    
    
#             # 2-1 단계: Input Molecule Properties
#             # =================================

#             try:
#                 prop_raw : QuerySet = io_mol_property(mol_id)
    
#             except PipelineError as e:
#                 _omit_mol.append(mol_id)
#                 print(f'Protein : {prot_id} - Molecule : {mol_id} : {e}')
#                 continue
            
#             # =================================
    
    
#             # 2-2 단계: Parse Fundamental Molecule Properties    
#             # =================================

#             try:
#                 prop : dict = parse_mol_property(prop_raw)
    
#             except PipelineError as e:
#                 _omit_mol.append(mol_id)
#                 print(f'Protein : {prot_id} - Molecule : {mol_id} : {e}')
#                 continue
    
#             prop_dict.update(prop)

#             # =================================
    
    
#             # 2-3 단계: Calculate Additional Properties
#             # =================================

#             canonical_smiles : str = prop_dict['canonical_smiles']
    
#             # 2-3-1 단계: Van Der Waals Volume
#             try:
#                 vdw : float = calc_vdwvolume(canonical_smiles)
    
#             except PipelineError as e:
#                 _omit_mol.append(mol_id)
#                 print(f'Protein : {prot_id} - Molecule : {mol_id} : {e}')
#                 continue
                
#             # 2-3-2 단계: Acidity, Basicity
#             # 미정

#             # =================================
            
#     OMIT_MOL[prot_id] = _omit_mol

#     # ============================
# ============================

#endregion


# test
if __name__ == '__main__':
    print('[0-0-0] Test: Get started')
    print()

    # Test
    # =========================================
    test_id = 'CHEMBL203'
    
    # level 1
    raw = io_prot_bindmol(test_id)
    molecule = raw[12837]
    print('[1-0-0] Test:')
    print(molecule)

    # level 2-0
    mol = parse_prot_query(molecule)                            # ID, dict(Properties) - Tuple
    print('[2-0-0] Test:')
    print(mol)

    # level 2-1
    mol_raw = io_mol_property(mol[0])                           # Molecule ID 해당하는 모든 Chemical Properties Input - Dict
    print('[2-1-0] Test:')
    print(mol_raw)

    # level 2-2
    mol_prop = parse_mol_property(mol_raw)                      # 필요한 Property만 모아서
    mol[1].update(mol_prop)                                     # Molecule Tuple의 Property Dict에 추가
    print('[2-2-0] Test:')
    print(mol)

    # level 2-3-1
    mol[1]['vdw'] = calc_vdwvolume(mol[1]['canonical_smiles'])  # Van der Waals Volume 계산 후 더해주기
    print('[2-3-1] Test:')
    print(mol)

    # level 2-3-2
    # =========================================
