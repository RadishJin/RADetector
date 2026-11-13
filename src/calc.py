#region MODULE & LIBRARY
# =========================================
from __future__ import annotations

# Standard Libraries

# Third-Party Libraries
from rdkit import Chem
from rdkit.Chem import AllChem, rdchem
from rdkit.Chem.SaltRemover import SaltRemover
from dimorphite_dl import DimorphiteDL as DD

# Custom Modules
from exceptions import EmbedFailError, NotValidError

# =========================================
#endregion
#region GLOBAL VARIABLE

_SALT_REMOVER = SaltRemover()
_PARAM = AllChem.ETKDGv3()   
_UNIQUE_IONIZER = DD()
_DIV_IONIZER = DD(min_ph = 3.9, max_ph = 9.1)

#endregion
#region INTERNAL

class _OptError(Exception):
    """ 최적화 실패 반환용 예외 """

def _parse_smiles(canonical_smiles : str) -> rdchem.Mol:
    """ smiles를 전용 class로 파싱 """

    mol = Chem.MolFromSmiles(canonical_smiles)

    if mol is None:
        raise NotValidError('잘못된 형태의 Canonical smiles')

    return mol
# Errors : [NotValidError]

def _deparse_mol(mol : rdchem.Mol) -> str:
    """ Mol 객체를 smiles로 전환 """

    if mol is None:
        raise NotValidError('None 타입의 Mol 객체가 입력됨')
    
    try: 
        smiles = Chem.MolToSmiles(mol)
    except RuntimeError:
        raise NotValidError('망가진 Mol 객체가 입력됨') 

    return smiles
# Errors : [NotValidError]

def _desalt_mol(mol : rdchem.Mol) -> rdchem.Mol:
    """ 염 제거 & Sanitize 함수 """

    try:
        strip_mol = _SALT_REMOVER.StripMol(mol, dontRemoveEverything= True)
    except RuntimeError:
        raise EmbedFailError('염 제거 실패')

    if strip_mol.GetNumAtoms() == 0:
        raise EmbedFailError('염 제거 후 분자 전체 소실')

    raw_n = mol.GetNumHeavyAtoms() 
    strip_n = strip_mol.GetNumHeavyAtoms()

    if raw_n > 5 and strip_n < (raw_n * 0.5):
        raise EmbedFailError('염 제거 후 중원자 반 이상 소실')

    try:
        Chem.SanitizeMol(strip_mol)
    except Chem.MolSanitizeException:
        raise EmbedFailError('염 제거 후 구조 재표준화 실패')

    return strip_mol
# Errors : [EmbedFailError]

def _opt_mol(mol : rdchem.Mol, maxiter : int, ise : bool = False) -> None | float:
    """ MMFF기반 물리화학 구조 최적화 """

    try:
        idx = AllChem.MMFFOptimizeMolecule(mol, maxIters = maxiter)        

        if idx != 0:
            raise _OptError
  
        if ise:
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant = "MMFF94")
            ff = AllChem.MMFFGetMoleculeForceField(mol, props)
            return ff.CalcEnergy()
        
    except (RuntimeError, _OptError):
        raise EmbedFailError('최적화 시도 실패')
# Errors : [EmbedFailError]

def _embed_mol_3d(mol : rdchem.Mol, maxseed : int, maxiter : int) -> rdchem.Mol:
    """ 3D 임베딩 """

    mol_std = Chem.AddHs(mol)
    success_list = []

    for seed in range(maxseed):
        psdmol = rdchem.Mol(mol_std)
        _PARAM.randomSeed = seed
        idx = AllChem.EmbedMolecule(psdmol, _PARAM)

        if idx != 0:
            continue

        try:
            success_list.append((_opt_mol(psdmol, maxiter, ise = True), psdmol))
        except EmbedFailError:
            continue

    if not success_list:
        raise EmbedFailError('최적화 성공 0/10 회')

    return min(success_list, key = lambda x: x[0])[1]
# Errors : [EmbedFailError]

#endregion
#region MAIN

def clean_smiles(raw_smiles : str) -> str:
    """ smiles의 염 제거, Sanitize하여 표준화하는 함수"""

    mol = _parse_smiles(raw_smiles)
    clean_mol = _desalt_mol(mol)
    clean_smiles = _deparse_mol(clean_mol)

    return clean_smiles

def calc_vdwvolume(canonical_smiles : str, maxseed : int, maxiter : int) -> float:
    """ smiles 를 이용하여 Van Der Waals Volume을 계산해주는 함수 """

    mol = _parse_smiles(canonical_smiles)

    mol = _embed_mol_3d(mol, maxseed, maxiter)

    return AllChem.ComputeMolVolume(mol)

def calc_formalcharge(canonical_smiles : str) -> float:
    """ pH 6.4에서의 Formal Charge를 계산하는 함수 """

    candidate_list = _UNIQUE_IONIZER.protonate(canonical_smiles)

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

    return sum(cand_charges) / len(cand_charges)

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