#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries
from dataclasses import dataclass, field

# Third-Party Libraries

# Custom Modules
from calc import clean_smiles, calc_vdwvolume

#endregion
#region MAIN

@dataclass
class Molecule:
    """ 단일 Molecule의 관련한 정보를 정리하는 클래스 """
    id : str
    smiles : str
    properties : dict = field(default_factory = dict)

    @classmethod
    def extended(cls,
                mol_id : str,
                raw_smiles : str,
                properties : dict = field(default_factory = dict),
                maxseed : int = 10,
                maxiter : int = 1000) -> Molecule:
        
        smiles = clean_smiles(raw_smiles)
        vdw : float = calc_vdwvolume(smiles, maxseed, maxiter)

        properties['vdw'] = vdw

        return cls(mol_id, smiles, properties)

#endregion
