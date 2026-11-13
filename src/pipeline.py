#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries

# Third-Party Libraries
from chembl_webresource_client.new_client import QuerySet

# Custom Modules & Classes
from io_parser import io_prot_bindmol, parse_prot_query, io_mol_property, parse_mol_property
from calc import clean_smiles, calc_vdwvolume
from config import PipelineError, Setting, Config
from domain import Molecule

#endregion
#region INTERNAL



#endregion
#region MAIN

setting = Setting()
config = Config(setting.target_proteins)

for prot_id in config.prot_ids:

    try:
        mol_raw : QuerySet = io_prot_bindmol(prot_id)

    except PipelineError as e:
        config.omits[prot_id] = 'All'
        print(f'Protein : {prot_id} : {e}')
        continue 

    config.omits[prot_id] = []
    for raw in mol_raw:
        data : tuple = parse_prot_query(raw)
        mol = Molecule(data[0], data[1], data[2])

        try:
            prop_raw : QuerySet = io_mol_property(mol.id)

        except PipelineError as e:
            config.omits[prot_id].append(mol.id)
            print(f'Protein : {prot_id} - Molecule : {mol.id} : {e}')
            continue

        try:
            prop : dict = parse_mol_property(prop_raw)

        except PipelineError as e:
            config.omits[prot_id].append(mol.id)
            print(f'Protein : {prot_id} - Molecule : {mol.id} : {e}')
            continue

        mol.properties.update(prop)

        try:
            vdw : float = calc_vdwvolume(mol.smiles)

        except PipelineError as e:
            config.omits[prot_id].append(mol.id)
            print(f'Protein : {prot_id} - Molecule : {mol.id} : {e}')
            continue

        try:
            fc = calc_fc(mol.smiles)

        except PipelineError as e:
            config.omits[prot_id].append(mol.id)
            print(f'Protein : {prot_id} - Molecule : {mol.id} : {e}')
            continue

        try:
            sigma = calc_sigma(mol.smiles)

        except PipelineError as e:
            config.omits[prot_id].append(mol.id)
            print(f'Protein : {prot_id} - Molecule : {mol.id} : {e}')
            continue

        mol.properties['vdw'] = vdw
        mol.properties['formal_charge'] = fc
        mol.properties['charge_sigma'] = sigma

        config.res_mols.append(mol)

        # 아래로는 미정

# 초안
def radetector():
    for prot_id in PROTEIN_IDS:
        _omit_mol : list = []

        # 1 : Input Binding Molecules
        try:
            mol_raw : QuerySet = io_prot_bindmol(prot_id)

        except PipelineError as e:
            OMIT_PROT.append(prot_id)
            print(f'Protein : {prot_id} : {e}')
            continue

        # 2-0 : Parse & Create Molecule Property Dictionary
        tot_mol_dict = {}
        for raw in mol_raw:
            mol : tuple = parse_prot_query(raw)
            mol_id : str = mol[0]
            prop_dict : dict[str, float] = mol[1]
    
            # 2-1 : Input Molecule Properties
            try:
                prop_raw : QuerySet = io_mol_property(mol_id)
            except PipelineError as e:
                _omit_mol.append(mol_id)
                print(f'Protein : {prot_id} - Molecule : {mol_id} : {e}')
                continue
    
            # 2-2 : Parse Fundamental Molecule Properties    
            try:
                prop : dict = parse_mol_property(prop_raw)
            except PipelineError as e:
                _omit_mol.append(mol_id)
                print(f'Protein : {prot_id} - Molecule : {mol_id} : {e}')
                continue
            prop_dict.update(prop)

            # 2-3-0 : Calculate Additional Properties
            prop_dict['canonical_smiles'] = clean_smiles(prop_dict['canonical_smiles'])
            canonical_smiles = prop_dict['canonical_smiles']
    
            # 2-3-1 : Van Der Waals Volume
            try:
                vdw : float = calc_vdwvolume(canonical_smiles)
    
            except PipelineError as e:
                _omit_mol.append(mol_id)
                print(f'Protein : {prot_id} - Molecule : {mol_id} : {e}')
                continue
                
            # 2-3-2 : fc, sigma
            fc = calc_fc(canonical_smiles)
            sigma = calc_sigma(canonical_smiles)

            # 2-4 : Final dict update
            prop_dict['vdw'] = vdw
            prop_dict['formal_charge'] = fc
            prop_dict['charge_sigma'] = sigma

            # 2-5 Upstream dict
            tot_mol_dict[mol_id] = prop_dict  

        OMIT_MOL[prot_id] = _omit_mol

        # 3-0 : Create multi-molecule array
        res_array = []#ndarray

        # 3-1 : Add each array into one array
        for mol_id, prop_dict in tot_mol_dict.items():
            res_array.append(array_mol(prop_dict))

        # 3-2 : transform to standard array
        array_standard(res_array)

        # 3-3 : (어떤 처리 할지 미정)
        array_transform(res_array)

        # 3-3 : Upstream dict
        RES_DICT[prot_id] = res_array

def prot_how_close(prot_id1, prot_id2):

    array1 = RES_DICT[prot_id1]
    array2 = RES_DICT[prot_id2]

    sim = how_similar(array1, array2)

    return sim





#endregion

# Test
# =========================================
test_id = 'CHEMBL203'

# 1
raw = io_prot_bindmol(test_id)
molecule = raw[12]   # 단일 Mol 조절
nnn = [raw[n] for n in range(100)]

# level 2-0
mol = parse_prot_query(molecule)                            # ID, dict(Properties) - Tuple
mol_list = []
for m in nnn:
    mol_list.append(parse_prot_query(m))

# level 2-1
mol_raw = io_mol_property(mol[0])                           # Molecule ID 해당하는 모든 Chemical Properties Input - Dict

# level 2-2
mol_prop = parse_mol_property(mol_raw)                      # 필요한 Property만 모아서
mol[1].update(mol_prop)                                     # Molecule Tuple의 Property Dict에 추가

# level 2-3-0
canonical_smiles = clean_smiles(mol[1]['canonical_smiles'])

# level 2-3-1
mol[1]['vdw'] = calc_vdwvolume(canonical_smiles)  # Van der Waals Volume 계산 후 더해주기
# print('[2-3-1] Test:')
print(mol)

# level 2-3-2
# =========================================
