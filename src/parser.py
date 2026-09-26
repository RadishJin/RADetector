import chembl_webresource_client.new_client as nc
import pandas as pd

# id_list = ["CHEMBL203", "CHEMBL1824", "CHEMBL1844", "CHEMBL301", "CHEMBL205"]

# Chembl data 라이브러리 이용해서 파싱 - 분자 구조 + 결합 정도

act = nc.new_client.activity
data = act.filter(
    target_chembl_id = "CHEMBL203",
    assay_type = 'B',
    standard_relation = "=",
    pchembl_value__isnull = False,
    canonical_smiles__isnull = False
).only([
    "molecule_chembl_id",
    "canonical_smiles",
    "pchembl_value"
])
df = pd.DataFrame(data)
print(df)

# SQLite 사용을 익혀야할듯
# pchembl value로 뭉개서 데이터를 받아오면서, Ki Kd IC50 같은 값의 출처 데이터 타입에 대한 정보가 누락됨


# test
# sample_keys = list(act.filter(target_chembl_id="CHEMBL203")[0].keys())
# print(sample_keys)

