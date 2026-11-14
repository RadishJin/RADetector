#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries

# Third-Party Libraries
import pandas as pd
import numpy as np
import ot

#endregion
#region GLOBAL VARIABLE

_FEATURE_COLS = [
    'vdw', 'sigma', 'aromatic_rings', 
    'hba', 'hbd', 'mw_freebase', 
    'psa', 'alogp', 'rtb'
]

#endregion
#region INTERNAL

def _zip_array(array : pd.DataFrame) -> tuple[np.ndarray]:
    """ 우리가 다룰 형태로 미리 전환 """
    coords = array[_FEATURE_COLS].values.astype(np.float64)
    raw_weights = array['pchembl'].values.astype(np.float64)
    weights = raw_weights / np.sum(raw_weights)

    return (coords, weights)

def _stdscaler(prot_a : pd.DataFrame, prot_b : pd.DataFrame) -> tuple[np.ndarray]:
    """ numpy 이용 두 좌표계 표준화 """
    a_coord, a_weight = _zip_array(prot_a)
    b_coord, b_weight = _zip_array(prot_b)

    comb = np.vstack([a_coord, b_coord])
    mean = np.mean(comb, axis=0)
    std = np.std(comb, axis=0)

    std[std == 0.0] = 1.0

    a_coord = (a_coord - mean) / std
    b_coord = (b_coord - mean) / std

    return a_coord, a_weight, b_coord, b_weight

#endregion
#region MAIN

def unit_array(prop : dict, pchembl : float, id : str) -> dict:
    return {
        'mol_id' : id,
        'pchembl' : pchembl,
        **prop  
    }

def build_array(unit_list : list[dict]) -> pd.DataFrame:
    return pd.DataFrame(unit_list)

def how_similar(prot_a : pd.DataFrame, prot_b : pd.DataFrame, gamma : float = 1.0) -> float:
    """ 두 단백질의 리간드 분포 유사도를 가장 가까우면 1 아예 다르면 0으로 1과 0 사이의 연속적인 값을 반환 """
    a_coord, a_weight, b_coord, b_weight = _stdscaler(prot_a, prot_b)

    cost_matrix = ot.dist(a_coord, b_coord, metric = 'euclidean')
    w_distance = ot.emd2(a_weight, b_weight, cost_matrix)
    sim = np.exp(-gamma * w_distance)

    return float(sim)

#endregion