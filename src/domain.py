#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries
from dataclasses import dataclass, field

# Third-Party Libraries

# Custom Modules
from calc import clean_smiles

#endregion
#region MAIN

@dataclass
class Molecule:
    """ 단일 Molecule의 관련한 정보를 정리하는 클래스 """
    id : str
    smiles : str
    properties : dict = field(default_factory = dict)

    def __post_init__(self):
        self.smiles = clean_smiles(self.smiles)

#endregion
