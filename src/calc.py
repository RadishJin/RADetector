from rdkit import Chem
from rdkit.Chem import AllChem, rdchem

from config import EmbedFailError, NotValidError

# ========================================= 내부 함수
# smiles를 전용 class로 파싱
def _parse_smiles(canonical_smiles : str) -> rdchem.Mol:

    mol = Chem.MolFromSmiles(canonical_smiles)  # RDkit 전용 Class로 Smiles 파싱

    if mol is None:
        raise NotValidError(f'{canonical_smiles} : 잘못된 형태의 Canonical smiles')

    return mol


# 3D 임베딩
def _embed_mol_3d(mol : rdchem.Mol) -> rdchem.Mol:

    mol_h = Chem.AddHs(mol)                     # 일괄적으로 Proton 붙이기
    
    param = AllChem.ETKDGv3()                   # 3D Embedding 전 ETKDGv3 기반 파라미터 설정
    param.randomSeed = 42                       # 재현성 확보 위한 시드 고정

    res = AllChem.EmbedMolecule(mol_h, param)   # 3D Embedding, 변수를 직접수정, 함수는 성공(0) 실패(-1) 반환
    if res < 0:
        raise EmbedFailError(f'3D 임베딩 실패')

    AllChem.UFFOptimizeMolecule(mol_h)          # 물리화학 기반 위치 Optimization, 변수 직접수정

    return mol_h
# =========================================


# =========================================
# Van Der Waals Volume calculation
def calc_vdwvolume(canonical_smiles : str) -> float:
    """ smiles 를 이용하여 Van Der Waals Volume을 계산해주는 함수 """

    mol = _parse_smiles(canonical_smiles)      # Smiles 클래스변환
    
    mol_h = _embed_mol_3d(mol)                 # 3D Embedding

    return AllChem.ComputeMolVolume(mol_h)
# =========================================


# =========================================
# pKa calculation
# dimorphite-DL로 라벨링 -> 베이스 애시드만 따로 서로 세기 비교 가능하게 계산
# =========================================
