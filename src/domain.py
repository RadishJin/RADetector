#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries
from dataclasses import dataclass, field
from typing import Any

# Custom Modules
from io_parser import io_mol_property_batch, parse_mol_property
from chem import clean_smiles, calc_vdwvolume, calc_sigmacharge
from build import unit_array, build_array, how_similar

#endregion
#region MAIN

def get_known_props(mol_ids : list, batch_size : int, max_retries : int) -> dict[str, dict]:
    """ 이미 알려져서 ChEMBL에 올라와있는 Molecule Property를 가져오는 함수 """
    prop_raw = io_mol_property_batch(mol_ids, batch_size, max_retries) 

    res = {}
    for mol_id in mol_ids:
        if mol_id in prop_raw and prop_raw[mol_id]:
            res[mol_id] = parse_mol_property(prop_raw[mol_id])

    return res

@dataclass
class Molecule:
    """ 단일 Molecule의 관련한 정보를 정리하는 클래스 """
    id : str
    smiles : str
    pchembl : float
    properties : dict = field(default_factory = dict)

    @classmethod
    def extended(cls,
                mol_id : str,
                raw_smiles : str,
                properties : dict = field(default_factory = dict),
                maxseed : int = 10,
                maxiter : int = 1000) -> Molecule:
        
        smiles : str = clean_smiles(raw_smiles)
        vdw : float = calc_vdwvolume(smiles, maxseed, maxiter)
        sigma : float = calc_sigmacharge(smiles)

        properties['vdw'] = vdw
        properties['sigma'] = sigma

        return cls(mol_id, smiles, properties)

@dataclass
class RADetect:
    """ Molecule 분포정보를 정리하는 클래스 """
    prot_id : str
    array : Any

    @classmethod
    def extended(cls,
                prot_id : str,
                mol_list : list[Molecule]) -> RADetect:

        pre = [unit_array(mol.properties, mol.pchembl, mol.id) for mol in mol_list]
        array = build_array(pre)

        return cls(
            prot_id = prot_id,
            array = array
        )

    def similar(self, target : RADetect, gamma : float = 1.0) -> float:
        return how_similar(self.array, target.array, gamma)

#endregion
