from io_parser import io_prot_bindmol, parse_prot_query, io_mol_property, parse_mol_property
from calc import calc_vdwvolume

from io_parser import NotFoundError, NotValidError
from calc import EmbedFailError

# 전역변수
PROTEIN_IDS = [
    'CHEMBL203',    # EFGR
    'CHEMBL1824',   # HER2
    'CHEMBL1844',   # CSF1R
    'CHEMBL301',    # CDK2
    'CHEMBL205'     # CA2
]

# test pipeline
test_id = 'CHEMBL203'

raw = io_prot_bindmol(test_id)          # Protein ID 해당하는 모든 Binding Molecules Input - QuerySet >> 여기서 NotFoundError 발생 가능성
                                                                                                 
molecule = raw[129]                     # 아무 Molecule이나 뽑아서

mol = parse_prot_query(molecule)        # ID, dict(Properties) - Tuple

mol_raw = io_mol_property(mol[0])       # Molecule ID 해당하는 모든 Chemical Properties Input - Dict  >> 여기서 NotFoundError 발생 가능성

mol_prop = parse_mol_property(mol_raw)  # 필요한 Property만 모아서                                     >> 여기서 NotValidError 발생 가능성

mol[1].update(mol_prop)                 # Molecule Tuple의 Property Dict에 추가

vdw = calc_vdwvolume(mol[1]['canonical_smiles'])                                                 # >> 여기서 EmbedFailError 발생 가능성
mol[1]['vdw'] = vdw                     # Van der Waals Volume 계산 후 더해주기

print(mol)

