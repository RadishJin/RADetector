########################## MODULE & LIBRARY ####
from __future__ import annotations

# Standard Libraries
# ============================
from typing import Iterable
# ============================

# Third-Party Libraries
# ============================
from chembl_webresource_client.new_client import new_client
# ============================

# Custom Modules
# ============================
from config import NotFoundError, NotValidError
# ============================
#### MODULE & LIBRARY ##########################


########################### GLOBAL VARIABLE ####
_TARGET_KEYS = {
        'alogp', 'psa', 'mw_freebase', 'rtb',
        'hba', 'hbd', 'aromatic_rings'
    }
# 필요한 Property Whitelist
#### GLOBAL VARIABLE ###########################


################################## INTERNAL ####
# ============================1
def _validate_query(query: Iterable[dict]) -> Iterable[dict]:
    """ 이터레이터(쿼리) 무결성 검증 함수. 값이 하나라도 있으면 True """

    try:
        _ = query[0]
# 쿼리 API 로딩 시도
    except (IndexError, TypeError):
        raise NotFoundError('잘못된 쿼리 요청 ID')
# [NotFoundError] : 쿼리에 담긴 값이 없을 시 예외처리

    return query
# 값이 하나라도 있는것이 확인된 쿼리 반환
# Errors : [NotFoundError]
# ============================1
#### INTERNAL ##################################


###################################### MAIN ####
# ============================1
def io_prot_bindmol(prot_id: str) -> Iterable[dict]:
    """ ChEMBL Binding Assay 데이터 가져오는 함수 """

    activity = new_client.activity
# Activity 층에 접근 가능한 객체 생성, 여기에 바인딩 어세이 존재
    query = activity.filter(
        target_chembl_id = prot_id,
        assay_type = 'B',
        pchembl_value__isnull = False
    ).only([
        'molecule_chembl_id',
        'canonical_smiles',
        'pchembl_value',
    ])
# SQL식 데이터 가져오기
    query = _validate_query(query)
# [NotFoundError] : 쿼리 무결성 검사, API 호출 1회 발생

    return query
# 무결한 쿼리 반환
# Errors : [NotFoundError]
# ============================1

# ============================2
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
# 딕셔너리 직접 가져와서 반환
# Errors : [NotFoundError]
# ============================2

# ============================3
# 수정 용이한 형태로 원데이터 파싱
def parse_prot_query(molecule: dict[str, str | float]) -> tuple[str, dict[str, str | float]]:
    """ ChEMBL Binding  Assay 데이터 정렬, 파싱 함수. 개별 데이터를 다룸. """

    return (molecule['molecule_chembl_id'], {k: v for k, v in molecule.items() if k != 'molecule_chembl_id' and k != 'value'})
# 뒤죽박죽 딕셔너리를 약속된 형태로 변환
# ============================3

# ============================4
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
# 원하는 property만 담긴 딕셔너리 반환
# Errors : [NotFoundError], [NotValidError]
# ============================4
#### MAIN ######################################