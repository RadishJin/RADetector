#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries
from dataclasses import dataclass, field
import logging
import os
from datetime import datetime

#endregion
#region MAIN

def setup_global_logger(name : str = "RADetector", log_dir : str = 'logs') -> logging.Logger:
    """
    전역 로거를 초기화하여 세팅함
    - Console : INFO level 이상의 log 출력
    - File : DEBUG 레벨 이상의 log 출력
    """
    os.makedirs(log_dir, exist_ok = True)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
        fmt="[%(asctime)s] [%(levelname)s] [%(name)s:%(lineno)d] - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)

    today_str = datetime.now().strftime("%Y%m%d")
    log_file_path = os.path.join(log_dir, f"{name.lower()}_{today_str}.log")
    
    file_handler = logging.FileHandler(log_file_path, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger

@dataclass
class Setting:
    """ 파이프라인 수정 가능 값들을 모아놓은 클래스 """

    # Molecule 3D Optimization 간 사용하는 시드는 0부터 순서대로 숫자를 쓰기로 했는데, 그 값의 상한선. 사실상 시드 개수
    maxseed : int = 10

    # Molecule 3D Optimization 간 최적 구조를 찾기 위해 연산을 몇 번 반복할지
    maxiter : int = 1000
    
    # Simularity Gamma 값
    gamma : float = 1.0

    # RADetect 산출할 Protein ID 목록
    target_proteins : list[str] = field(
        default_factory = lambda: [
        'CHEMBL203',    # EFGR
        'CHEMBL1824',   # HER2
        'CHEMBL1844',   # CSF1R
        'CHEMBL301',    # CDK2
        'CHEMBL205'     # CA2
        ]
    )

    # API에 한번에 요청할 Batch Size
    batch_size : int = 50

    # API 요청 재시도 횟수
    max_retries : int = 5

@dataclass
class Config:
    """ 데이터 총괄 흐름 관리용 클래스 """
    omits : dict = field(default_factory = dict)
    res_prots : dict = field(default_factory = dict)
    res_mols : list = field(default_factory = list)

#endregion