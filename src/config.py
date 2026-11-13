#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries
from dataclasses import dataclass, field

#endregion
#region MAIN

@dataclass
class Setting:
    """ 파이프라인 수정 가능 값들을 모아놓은 클래스 """

    # Molecule 3D Optimization 간 사용하는 시드는 0부터 순서대로 숫자를 쓰기로 했는데, 그 값의 상한선. 사실상 시드 개수
    maxseed : int = 10

    # Molecule 3D Optimization 간 최적 구조를 찾기 위해 연산을 몇 번 반복할지
    maxiter : int = 1000

    # Molecule Property Parsing 간 Whitelist
    target_properties : set = field(
        default_factory = lambda: {
        'alogp', 'psa', 'mw_freebase', 'rtb',
        'hba', 'hbd', 'aromatic_rings'
        }
    )

    # RADetect 산출할 Protein ID 목록
    target_proteins : list[str] = field(
        default_factory = lambda: [
        'CHEMBL203',    # EFGR
        'CHEMBL1824',   # HER2
        'CHEMBL1844',   # CSF1R
        'CHEMBL301',    # CDK2
        'CHEMBL205'     # CA2
        ]
    )

@dataclass
class Config:
    """ 데이터 총괄 흐름 관리용 클래스 """
    prot_ids : list[str]
    omits : dict = field(default_factory = dict)
    res_prots : dict = field(default_factory = dict)
    res_mols : list = field(default_factory = list)

#endregion