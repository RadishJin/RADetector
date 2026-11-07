from __future__ import annotations

from chembl_webresource_client.new_client import new_client
from typing import Iterable

from config import NotFoundError, NotValidError


# 내부 함수
# =========================================
# 이터레이터(쿼리) 무결성 검증 함수. 값이 하나라도 있으면 True
def _validate_query(query: Iterable[dict]) -> Iterable[dict]:

    try:
        _ = query[0]
    except (IndexError, TypeError):
        raise NotFoundError('잘못된 쿼리 요청 ID')

    return query
# =========================================


# =========================================
def io_prot_bindmol(prot_id: str) -> Iterable[dict]:
    """ ChEMBL Binding Assay 데이터 가져옴 """

    # Activity 층 접근 -> 여기에 바인딩 어세이 존재
    activity = new_client.activity

    # SQL식 데이터 가져오기
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


# Chemical Property 
def io_mol_property(mol_id: str) -> dict:
    """ ChEMBL Molecule Property 데이터 가져옴 """

    # Molecule 층 접근 -> 여기에 pKa 존재
    molecule = new_client.molecule

    # SQL식 데이터 가져오기
    query = molecule.filter(
        molecule_chembl_id = mol_id
    ).only([
        'molecule_properties',
    ])

    query = _validate_query(query)

    return query[0]['molecule_properties']



# 수정 용이한 형태로 원데이터 파싱
def parse_prot_query(molecule: dict[str, str | float]) -> tuple[str, dict[str, str | float]]:
    """ ChEMBL Binding  Assay 데이터 정렬, 파싱 함수. 개별 데이터를 다룸. """

    return (molecule['molecule_chembl_id'], {k: v for k, v in molecule.items() if k != 'molecule_chembl_id' and k != 'value'})


# property 데이터 파싱
def parse_mol_property(mol_prop: dict[str, str | float]) -> dict[str, float]:
    """ ChEMBL Molecule pKa 데이터 정렬, 파싱 함수. 개별 데이터를 다룸. """

    target_keys = {
        'alogp', 'psa', 'mw_freebase', 'rtb',
        'hba', 'hbd', 'aromatic_rings'
    }
    
    parsed = {}
    for k in target_keys:
        val = mol_prop.get(k)

        if val is None:
            raise NotFoundError(f'{k}값 없음')
        
        try:
            parsed[k] = float(val)
        except (ValueError, TypeError):
            raise NotValidError('잘못된 형태의 Property 값')
            
    return parsed
# =========================================