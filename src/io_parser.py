#region MODULE & LIBRARY
# =========================================
from __future__ import annotations

# Standard Libraries

# Third-Party Libraries
from chembl_webresource_client.new_client import new_client, QuerySet

# Custom Modules
from config import NotFoundError, NotValidError

# =========================================
#endregion
#region GLOBAL VARIABLE
# =========================================

_TARGET_KEYS = {
        'alogp', 'psa', 'mw_freebase', 'rtb',
        'hba', 'hbd', 'aromatic_rings'
    }

# =========================================
#endregion
#region INTERNAL
# =========================================

def _validate_query(query: QuerySet) -> QuerySet:
    """ 이터레이터(쿼리) 무결성 검증 함수. 값이 하나라도 있으면 True """

    try:
        _ = query[0]
    except (IndexError, TypeError):
        raise NotFoundError('잘못된 쿼리 요청 ID')

    return query
# Errors : [NotFoundError]

# =========================================
#endregion
#region MAIN
# =========================================

def io_prot_bindmol(prot_id: str) -> QuerySet:
    """ ChEMBL Binding Assay 데이터 가져오는 함수 """

    activity = new_client.activity

    query = activity.filter(
        target_chembl_id = prot_id,
        assay_type = 'B',
        pchembl_value__isnull = False
    ).only([
        'molecule_chembl_id',
        'canonical_smiles',
        'pchembl_value',
    ])

    query = _validate_query(query)

    return query

def io_mol_property(mol_id: str) -> dict:
    """ ChEMBL Molecule Property 데이터 가져옴 """

    molecule = new_client.molecule
# Molecule 층에 접근 가능한 객체 생성, 여기에 Properties 존재
    query = molecule.filter(
        molecule_chembl_id = mol_id
    ).only([
        'molecule_properties',
    ])
# SQL식 데이터 가져오기
    query = _validate_query(query)
# [NotFoundError] : 쿼리 무결성 검사, API 호출 1회 발생

    return query[0]['molecule_properties']

def parse_prot_query(molecule: dict[str, str | float]) -> tuple[str, str, dict[str, str | float]]:
    """ ChEMBL Binding  Assay 데이터 정렬, 파싱 함수. 개별 데이터를 다룸. """

    return (molecule['molecule_chembl_id'],
            molecule['canonical_smiles'],
            {k: v for k, v in molecule.items() if k != 'molecule_chembl_id' and k != 'value' and k != 'canonical_smiles'}
    )

def parse_mol_property(mol_prop: dict[str, str | float]) -> dict[str, float]:
    """ ChEMBL Molecule pKa 데이터 정렬, 파싱 함수. 개별 데이터를 다룸. """

    parsed = {}
    for k in _TARGET_KEYS:
        val = mol_prop.get(k)
# target property 입력 딕셔너리에서 불러오기 시도
        if val is None:
            raise NotFoundError(f'{k}값 없음')
# [NotFoundError] : 타겟 키가 존재하지 않는 상황 예외처리
        try:
            parsed[k] = float(val)
# 모든 값을 약속된 자료형(float)으로 변환하여 새로운 딕셔너리에 저장
        except (ValueError, TypeError):
            raise NotValidError('잘못된 형태의 Property 값')
# [NotValidError] : 약속된 자료형이 아닌 경우 예외처리
            
    return parsed

# =========================================
#endregion