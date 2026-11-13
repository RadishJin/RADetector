#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries

# Third-Party Libraries

# Custom Modules & Classes
from io_parser import io_mol_property, parse_mol_property

#endregion
#region MAIN

def get_known_props(mol_id : str) -> dict[str, float]:
    """ 이미 알려져서 ChEMBL에 올라와있는 Molecule Property를 가져오는 함수 """
    prop_raw = io_mol_property(mol_id)
    prop : dict = parse_mol_property(prop_raw)

    return prop

#endregion

