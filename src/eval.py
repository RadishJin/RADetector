import torch

# .pt 파일 불러오기, 변수에 저장까지
id_list = ["1CRN", "1CLL", "5DK3"]
data = {}
for id in id_list:
    ca = torch.load(f"test_dataset/{id}_ca.pt")
    bb = torch.load(f"test_dataset/{id}_bb.pt")
    data[f"{id}_ca"] = ca
    data[f"{id}_bb"] = bb

# print(data.keys())
# print(data.values())