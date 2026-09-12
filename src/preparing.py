from parsing import parsing                 # 파싱 데이터 생성
import biotite.structure.io.pdbx as pdbx    # 읽어오기      

# parsing()   # 파싱 데이터 생성. 이미 있으면 덮어써버림

id_list = ["1CRN", "1CLL", "5DK3"]
data_dict = {}
for id in id_list:    
    with open(f"raw_dataset/parsed_{id}.cif", "rt", encoding = "utf-8") as f:
        data_dict[f"{id}"] = pdbx.CIFFile.read(f)

# test
# print(data_dict["1CRN"])


# 각 target의 특징에 따라 torsion angle 기반 decoy 생성

# 1. 1CRN Crambin (가장 기본형)
# (1) Gaussian Coordiante Noise

# (2) Global Perturbation

# (3) Local Extreme Perturbation









