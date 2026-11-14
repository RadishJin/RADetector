#region MODULE & LIBRARY
from __future__ import annotations

# Custom Modules & Classes
from io_parser import io_prot_bindmol, parse_mol_query
from domain import get_known_props, Molecule, RADetect

from config import Config, Setting
from exceptions import PipelineError

#endregion
#region MAIN

def radetector() -> None:
    """ 미정 """
    

def how_similar() -> None:
    """ 미정 """

def pipeline(config : Config, setting : Setting) -> None:
    for prot_id in config.prot_ids:

        try:
            mol_raw = io_prot_bindmol(prot_id)
        except PipelineError as e:
            config.omits[prot_id] = 'All'
            print(f'Protein : {prot_id} : {e}')
            continue 

        config.omits[prot_id] = []
        for raw in mol_raw:
            data : tuple = parse_mol_query(raw)

            try:
                mol = Molecule.extended(data[0], data[1], data[2], data[3], setting.maxseed, setting.maxiter)

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

        radetect = MolArray.extended(prot_id, config.res_mols)
        
            # 아래로는 미정



#endregion