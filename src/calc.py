#region MODULE & LIBRARY

from __future__ import annotations

# Standard Libraries
# ============================
# ============================

# Third-Party Libraries
# ============================
from rdkit import Chem
from rdkit.Chem import AllChem, rdchem
from rdkit.Chem.SaltRemover import SaltRemover

from dimorphite_dl import DimorphiteDL as DD
# ============================

# Custom Modules
# ============================
from config import EmbedFailError, NotValidError
# ============================

#endregion
#region GLOBAL VARIABLE

_SALT_REMOVER = SaltRemover()
# 염 제거 전용 객체 생성
_PARAM = AllChem.ETKDGv3()
# 3D Embedding 전 ETKDGv3 기반 파라미터 설정      
_N_SEED = 10
# 시드로 이용할 수의 개수
_MAXITER = 1000
# 최적화 반복 계산 횟수
_UNIQUE_IONIZER = DD()
# 포멀차지 계산용 프로토네이션
_DIV_IONIZER = DD(min_ph = 3.9, max_ph = 9.1)
# 이온화 양상 계산용 프로토네이션

#endregion
#region INTERNAL

class _OptError(Exception):
    """ 최적화 실패 반환용 예외 """

def _parse_smiles(canonical_smiles : str) -> rdchem.Mol:
    """ smiles를 전용 class로 파싱 """

    mol = Chem.MolFromSmiles(canonical_smiles)
# rdchem.Mol Class로 Smiles 파싱 시도
    if mol is None:
        raise NotValidError('잘못된 형태의 Canonical smiles')
# [NotValidError] : 파싱 실패 예외처리
    return mol
# Errors : [NotValidError]

def _desalt_mol(mol : rdchem.Mol) -> rdchem.Mol:
    """ 염 제거 & Sanitize 함수 """

    try:
        strip_mol = _SALT_REMOVER.StripMol(mol, dontRemoveEverything= True)
# 원본 Mol에서 Salt 제거 시도
    except RuntimeError:
        raise EmbedFailError('염 제거 실패')
# [EmbedFailError] : 염 제거를 실패하는 경우 예외처리
    if strip_mol.GetNumAtoms() == 0:
        raise EmbedFailError('염 제거 후 분자 전체 소실')
# [EmbedFailError] : 염 제거 후 분자 전체 사라진 경우 예외처리
    raw_n = mol.GetNumHeavyAtoms() 
    strip_n = strip_mol.GetNumHeavyAtoms()
# 중원자(H 제외 전부) 개수를 염 제거 전후로 Count
    if raw_n > 5 and strip_n < (raw_n * 0.5):
        raise EmbedFailError('염 제거 후 중원자 반 이상 소실')
# [EmbedFailError] : 염 제거 후 중원자가 절반 이상 사라진 경우 예외처리
    try:
        Chem.SanitizeMol(strip_mol)
# 염 제거 이후 망가진 구조를 재표준화 시도
    except Chem.MolSanitizeException:
        raise EmbedFailError('염 제거 후 구조 재표준화 실패')
# [EmbedFailError] : 염 제거 이후 재표준화 실패한 경우 예외처리

    return strip_mol
# Errors : [EmbedFailError]

def _opt_mol(mol : rdchem.Mol, ise = False) -> None | float:
    """ MMFF기반 물리화학 구조 최적화 """

    try:
        idx = AllChem.MMFFOptimizeMolecule(mol, maxIters = _MAXITER)        
# MMFF 기반 위치 최적화(물리화학) 시도. 함수는 성공(0) 실패(1) 반환
        if idx != 0:
            raise _OptError
# MMFF가 최적해를 찾지 못할 시 [_OptError]  
        if ise:
# 사용자가 에너지 수치 반환을 요구한 경우 - 반환값 생김(에너지 상태)
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant = "MMFF94")
            ff = AllChem.MMFFGetMoleculeForceField(mol, props)
# MMFF 기반 Target Molecule의 ForceField 생성 후
            return ff.CalcEnergy()
# 에너지 계산하여 반환
    except (RuntimeError, _OptError):
        raise EmbedFailError('최적화 시도 실패')
# [EmbedFailError] : 최적화 시도 실패[RuntimeError], 최적화 중 최적해 계산 실패[_OptError] 모두 같은 예외처리
    
    return None
# Errors : [EmbedFailError]

def _embed_mol_3d(mol : rdchem.Mol) -> rdchem.Mol:
    """ 3D 임베딩 """

    mol_std = _desalt_mol(mol)
# [EmbedFailError] : 원본 분자 염 제거 시도
    mol_std = Chem.AddHs(mol_std)
# 암시적 수소를 명시적 수소로 생성
    success_list = []
    for seed in range(_N_SEED):
# 시드(0~9) 생성
        psdmol = rdchem.Mol(mol_std)
        _PARAM.randomSeed = seed
# 새로운 변수에 분자 클로닝, 시드 입력
        idx = AllChem.EmbedMolecule(psdmol, _PARAM)
# 3D Embedding 시도, 변수를 직접수정, 함수는 성공(0) 실패(-1) 반환
        if idx != 0:
            continue
# 성공이 아닌 모든 경우 폐기
        try:
            success_list.append((_opt_mol(psdmol, ise = True), psdmol))
# 임베딩 결과 최적화하여 에너지와 수정된 객체를 성공 리스트에 추가, 실패 시 [EmbedFailError]
        except EmbedFailError:
            continue
# 최적화 실패 시 폐기
    if not success_list:
        raise EmbedFailError('최적화 성공 0/10 회')
# [EmbedFailError] : 전량 최적화 실패 예외처리

    return min(success_list, key = lambda x: x[0])[1]
# Errors : [EmbedFailError]

#endregion
#region MAIN

def calc_vdwvolume(canonical_smiles : str) -> float:
    """ smiles 를 이용하여 Van Der Waals Volume을 계산해주는 함수 """

    mol = _parse_smiles(canonical_smiles)
# [NotValidError] : smiles 파싱
    mol = _embed_mol_3d(mol)
# [EmbedFailError] : 3D embedding

    return AllChem.ComputeMolVolume(mol)

def calc_formalcharge(canonical_smiles : str) -> float:
    """ pH 6.4에서의 Formal Charge를 계산하는 함수 """

    pre_mol = _parse_smiles(canonical_smiles)
    cand_mol = _desalt_mol(pre_mol)
    clean_smiles = Chem.MolToSmiles(cand_mol)
    candidate_list = _UNIQUE_IONIZER.protonate(clean_smiles)
    if not candidate_list:
        raise EmbedFailError('Protonation 실패')
    cand_charges = []
    for cand in candidate_list:
        try:
            cand_mol = _parse_smiles(cand)
        except NotValidError:
            continue
        cand_charges.append(sum([atom.GetFormalCharge() for atom in cand_mol.GetAtoms()]))
    if not cand_charges:
        raise EmbedFailError('유효한 Protonation 0건')

    print(len(cand_charges))

    return sum(cand_charges) / len(cand_charges)



    print(candidate_list)

def calc_sigmacharge(canonical_smiles : str) -> float:
    pass

#endregion

# test 
if __name__ == '__main__':
    from pipeline import mol_list

    for mol in mol_list:
        smiles = mol[1]['canonical_smiles']
        print(smiles)
        print(calc_formalcharge(smiles))

# 산염기는 어디로 가야하나
# 적정 시뮬레이터 전환: 단일 pH 상수를 버리고, pH 구간별 프로토머 앙상블을 스캔하는 적정 곡선 시뮬레이터 구축.
# Coarse-to-Fine 탐색: 넓은 범위를 대충 훑은 뒤, 전하가 급변하거나 0이 되는 변곡점(Inflection Point) 및 중성 도달 지점(0-Crossing)을 촘촘하게 타겟팅.
# 고유 임계 좌표 추출: 분자마다 전하가 꺾이는 고유의 임계 pH 좌표를 추출해 피처로 활용 (데이터 뭉개짐 원천 차단).