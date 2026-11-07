from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem, rdchem
from rdkit.Chem.SaltRemover import SaltRemover

from config import EmbedFailError, NotValidError


# 내부 변수
# =========================================
_salt_remover = SaltRemover()
_param = AllChem.ETKDGv3()                     # 3D Embedding 전 ETKDGv3 기반 파라미터 설정
_n_seed = 10
_maxiter = 1000
# =========================================


# 내부 함수&클래스
# =========================================
class _OptError(Exception): pass

# smiles를 전용 class로 파싱
def _parse_smiles(canonical_smiles : str) -> rdchem.Mol:

    mol = Chem.MolFromSmiles(canonical_smiles)  # RDkit 전용 Class로 Smiles 파싱

    if mol is None:
        raise NotValidError('잘못된 형태의 Canonical smiles')

    return mol

# 염 제거 & Sanitize 함수
def _desalt_mol(mol : rdchem.Mol) -> rdchem.Mol:

    try:
        strip_mol = _salt_remover.StripMol(mol, dontRemoveEverything= True)
    except RuntimeError:
        raise EmbedFailError('염 제거 중 끔찍한 오류')

    # 1차 정합성 검사 - 존재하는가?
    if strip_mol.GetNumAtoms() == 0:
        raise EmbedFailError('염 제거 후 분자 소실')

    # 2차 정합성 검사 - 반 이상 소실되었는가?
    raw_n = mol.GetNumHeavyAtoms()
    strip_n = strip_mol.GetNumHeavyAtoms()
    if raw_n > 5 and strip_n < (raw_n * 0.5):
        raise EmbedFailError('염 제거했더니 Thanos')

    # 3차 정합성 검사 - 다시 정상적인 구조로 돌아가는가?
    try:
        Chem.SanitizeMol(strip_mol)
    except Chem.MolSanitizeException:
        raise EmbedFailError('Sanitize 실패')

    return strip_mol      

# MMFF기반 물리화학 구조 최적화
def _opt_mol(mol : rdchem.Mol, ise = False) -> None | float:

    try:
        idx = AllChem.MMFFOptimizeMolecule(mol, maxIters = _maxiter)        # MMFF 기반 위치 최적화(물리화학). 함수는 성공(0) 실패(1) 반환

        if idx != 0:
            raise _OptError('MMFF 최적화 - 뒤틀린 구조')
        
        if ise:
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant = "MMFF94")
            ff = AllChem.MMFFGetMoleculeForceField(mol, props)

            return ff.CalcEnergy()
        
    except (RuntimeError, _OptError):
        raise EmbedFailError('MMFF 최적화 - 끔찍한 에러')
    
    return None

# 3D 임베딩
def _embed_mol_3d(mol : rdchem.Mol) -> rdchem.Mol:

    mol_std = _desalt_mol(mol)                     # 염 제거
    mol_std = Chem.AddHs(mol_std)                  # 일괄적으로 Proton 붙이기

    success_list = []
    for seed in range(_n_seed):
        psdmol = rdchem.Mol(mol_std)
        _param.randomSeed = seed

        idx = AllChem.EmbedMolecule(psdmol, _param)   # 3D Embedding, 변수를 직접수정, 함수는 성공(0) 실패(-1) 반환
        if idx != 0:
            continue

        try:
            success_list.append((_opt_mol(psdmol, ise = True), psdmol))
        except EmbedFailError:
            continue

    if not success_list:
        raise EmbedFailError('최적화 성공 0/10 회')

    return min(success_list, key = lambda x: x[0])[1]
# =========================================


# =========================================
def calc_vdwvolume(canonical_smiles : str) -> float:
    """ smiles 를 이용하여 Van Der Waals Volume을 계산해주는 함수 """

    mol = _parse_smiles(canonical_smiles)      # Smiles 클래스변환
    
    mol = _embed_mol_3d(mol)                 # 3D Embedding

    return AllChem.ComputeMolVolume(mol)
# =========================================


# =========================================
# pKa calculation
# dimorphite-DL로 라벨링 -> 베이스 애시드만 따로 서로 세기 비교 가능하게 계산
# =========================================


# test 
# from pipeline import mol
