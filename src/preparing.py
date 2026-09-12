from parsing import parsing                 # 파싱 데이터 생성
import biotite.structure.io.pdbx as pdbx    # 읽어오기     
import biotite.structure as struc
import torch                                # 텐서 변환, 노이즈 생성용 
from math import pi


# parsing()   # 파싱 데이터 생성. 이미 있으면 덮어써버림

id_list = ["1CRN", "1CLL", "5DK3"]
data_dict = {}
for id in id_list:    
    with open(f"raw_dataset/parsed_{id}.cif", "rt", encoding = "utf-8") as f:
        data_dict[f"{id}"] = pdbx.get_structure(pdbx.CIFFile.read(f))

# test
# print(data_dict["1CLL"])



# 각 target의 특징에 따라 torsion angle 기반 decoy 생성


# 1. 1CRN Crambin (가장 기본형)


# (1) Gaussian Coordiante Noise (sigma = 0.5, 1.0, 2.0, 5.0 angstrom)

bb_coord = struc.coord(data_dict["1CRN"])
bb_tensor = torch.tensor(bb_coord)
sigmas = [0.5, 1.0, 2.0, 5.0]
decoy_dict = {}
for sigma in sigmas:
    noise = torch.randn_like(bb_tensor) * sigma
    noise_coord = bb_tensor + noise
    noise_coord = noise_coord.detach().cpu().numpy()
    pre1 = data_dict["1CRN"].copy()
    struc.coord(pre1)[:] = noise_coord
    decoy_dict[f"1CRN_{sigma}angstrom"] = pre1
    
    
# test
# print(noise_coord)
# print(decoy_dict["1CRN_5.0angstrom"])
# print(decoy_dict.keys())


# (2) Global Perturbation (Gaussian, sigma = 5, 10 degree)

# 바꿀 토션앵글 값 구하기
bb_phi, bb_psi, bb_omega = struc.dihedral_backbone(data_dict["1CRN"])
bb_phi_tensor = torch.tensor(bb_phi)
bb_psi_tensor = torch.tensor(bb_psi)
bb_omega_tensor = torch.tensor(bb_omega)
tensor_list = [bb_phi_tensor, bb_psi_tensor, bb_omega_tensor]
angle_noise_list = [pi/36, pi/18]
pre_decoy = []
for angle in angle_noise_list:
    for tensor in tensor_list:
        noise = torch.randn_like(tensor) * angle
        noise = noise.detach().cpu().numpy().tolist()
        pre_decoy.append(noise[0])

# print(bb_omega)
# print(pre_decoy[5])

# 델타 토션앵글 리스트화
five_degree = [i for j in zip(pre_decoy[0], pre_decoy[1], pre_decoy[2]) for i in j]
ten_degree = [i for j in zip(pre_decoy[3], pre_decoy[4], pre_decoy[5]) for i in j]
# print(ten_degree)

# 처음과 끝 nan으로 바꾸기
five_degree = five_degree[1:-1]
ten_degree = ten_degree[1:-1]
# print(five_degree, ten_degree)

pre2 = data_dict["1CRN"][0].copy()
# print(pre[0].coord)

# N-term부터 하나씩 돌리면서 새 구조 만들기
k = 0
for angle in five_degree, ten_degree:
    for i in range(len(angle)):
        if i+2 > len(pre2):
            break
        axis = pre2[i+1].coord - pre2[i].coord
        support = pre2[i+1].coord
        downstream = pre2[i+2:]
        downstream = struc.rotate_about_axis(
            downstream,
            angle = angle[i],
            axis = axis,
            support = support
        )
        struc.coord(pre2[i+2:])[:] = struc.coord(downstream)
    k += 5
    decoy_dict[f"1CRN_{k}degree"] = pre2


# test
# print(data_dict["1CRN"][0])
# print(pre1[0][-5:])
# print(five_degree)
# print(pre2)

print(decoy_dict.keys())






# (3) Local Extreme Perturbation









