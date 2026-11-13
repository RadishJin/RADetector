#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries

# Third-Party Libraries
from chembl_webresource_client.new_client import QuerySet

# Custom Modules & Classes
from io_parser import io_prot_bindmol, parse_mol_query
from pipeline import get_known_props
from exceptions import PipelineError
from config import Setting, Config
from domain import Molecule

#endregion
#region MAIN

setting = Setting(target_proteins= ['CHEMBL203'])
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
        data : tuple = parse_mol_query(raw)

        try:
            mol = Molecule.extended(data[0], data[1], data[2], setting.maxseed, setting.maxiter)

        except PipelineError as e:
                    config.omits[prot_id].append(data[0])
                    print(f'Protein : {prot_id} - Molecule : {data[0]} : {e}')
                    continue
        
        try:
            mol.properties.update(get_known_props(mol.id))

        except PipelineError as e:
            config.omits[prot_id].append(mol.id)
            print(f'Protein : {prot_id} - Molecule : {mol.id} : {e}')
            continue

        print(f'Successful with {prot_id} - {mol.id}')
        config.res_mols.append(mol)
    
        # 아래로는 미정

print(len(mol_raw))
print(len(config.omits))
print(len(config.res_mols))

#endregion

# batch API
# MultiProcessing