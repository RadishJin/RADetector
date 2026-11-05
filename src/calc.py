from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem


# =========================================
# Custom Exception
class EmbedFailError(Exception):
    """ 분자 3D 임베딩 실패 """
# =========================================


# =========================================
# Van Der Waals Volume calculation
def calc_vdwvolume(canonical_smiles : Any) -> float:
    """ smiles 를 이용하여 Van Der Waals Volume을 계산해주는 함수 """

    mol = Chem.MolFromSmiles(canonical_smiles)  # RDkit 전용 Class로 Smiles 파싱
    if mol is None:
        raise EmbedFailError(f'{canonical_smiles} : Canonical smiles 파싱 실패')
    
    mol_h = Chem.AddHs(mol)                     # Volume 계산 전 일괄적으로 Proton 붙이기

    param = AllChem.ETKDGv3()                   # 3D Embedding 전 ETKDGv3 기반 파라미터 설정
    param.randomSeed = 42                       # 재현성 확보 위한 시드 고정

    res = AllChem.EmbedMolecule(mol_h, param)   # 3D Embedding, 변수를 직접수정, 함수는 성공(0) 실패(-1) 반환
    if res < 0:
        raise EmbedFailError(f'{canonical_smiles} : 3D 임베딩 실패')

    AllChem.UFFOptimizeMolecule(mol_h)          # 물리화학 기반 위치 Optimization, 변수 직접수정

    volume = AllChem.ComputeMolVolume(mol_h)    # Volume Calculate

    return volume
# =========================================