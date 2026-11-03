from chembl_webresource_client.new_client import new_client

class IDNotFoundError(ValueError):
    """ ChEMBL ID에 해당하는 바인딩 어세이 데이터 미존재 """

def io_chembl_b_molecule(prt_id: str) -> list[dict[str, str | float]]:
    """ ChEMBL Binding Assay 데이터 입력 함수 """

    # Activity 층 접근 -> 여기에 바인딩 어세이 존재
    activity = new_client.activity

    # SQL식 데이터 가져오기
    data = activity.filter(
        target_chembl_id = prt_id,
        assay_type = 'B',
        pchembl_value__isnull = False
    ).only([
        'molecule_chembl_id',
        'canonical_smiles',
        'pchembl_value',
    ])

    # 비어있는지 여부 예외처리
    data = list(data)
    if not data:
        raise IDNotFoundError(f'ChEMBL ID {prt_id}에 해당하는 바인딩 어세이 데이터가 존재하지 않습니다.')
    else:
        return data

def parse_chembl_b_molecule(raw: list[dict[str, str | float]]) -> dict[str, dict[str, str | float]]:
    """ ChEMBL Binding Assay 데이터 정렬, 파싱 함수 """

    return {m['molecule_chembl_id']: {k: v for k, v in m.items() if k != 'molecule_chembl_id'} for m in raw}

prt_id = ['CHEMBL203']
m_raw = io_chembl_b_molecule(prt_id)
m_dict = parse_chembl_b_molecule(m_raw)

#  재정렬




# test
# if __name__ == '__main__':
#     prt_id = ['CHEMBL203']
#     print(len(io_chembl_b_molecule(prt_id)))
