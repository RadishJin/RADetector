from chembl_webresource_client.new_client import new_client
from typing import Iterator


# 이터레이터 형태로 Input data 가져오는 함수
def io_chembl_b_molecule(prt_id: str) -> Iterator[dict]:
    """ ChEMBL Binding Assay 데이터를 이터레이터 형태로 가져옴 """

    # Activity 층 접근 -> 여기에 바인딩 어세이 존재
    activity = new_client.activity

    # SQL식 데이터 가져오기
    query = activity.filter(
        target_chembl_id = prt_id,
        assay_type = 'B',
        pchembl_value__isnull = False
    ).only([
        'molecule_chembl_id',
        'canonical_smiles',
        'pchembl_value',
    ])

    return query


# 수정 용이한 형태로 원데이터 파싱
def parse_chembl_b_molecule(raw: dict[str, str | float]) -> dict[str, dict[str, str | float]]:
    """ ChEMBL Binding Assay 데이터 정렬, 파싱 함수. 개별 데이터를 다룸. """

    return {raw['molecule_chembl_id']: {k: v for k, v in raw.items() if k != 'molecule_chembl_id' and k != 'value'}}


# test
if __name__ == '__main__':
    prt_id = 'CHEMBL203'
    m_raw = io_chembl_b_molecule(prt_id)
    m_test = parse_chembl_b_molecule(next(m_raw))
    print(m_test)