#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries
import logging

# Custom Modules & Classes
from pipeline import radetector
from domain import RADetect
from config import Setting, Config, setup_global_logger
from exceptions import PipelineError

#endregion
#region GLOBAL VARIABLE

logger = setup_global_logger()
setting = Setting()
config = Config()

#endregion
#region MAIN

# ChEMBL 클라이언트나 urllib3 쪽 로그는 너무 시끄러우니 WARNING 이상만 찍도록 잠재우기
logging.getLogger("chembl_webresource_client").setLevel(logging.WARNING)

def main():
    logger.info("===== RADetector 가동 시작 =====")
    rad_list = []
    for prot_id in setting.target_proteins:
        config.omits[prot_id] = []
        
        try:
            rad_list.append(radetector(prot_id, setting, config))
            logger.info(f'===== Protein : {prot_id} - RADetect 생성 성공 =====')

        except PipelineError as e:
            config.omits[prot_id] = 'All'
            logger.error(f'===== 실패한 Protein : {prot_id} - {e} =====')
            continue 

    if not (len(rad_list) == len(setting.target_proteins)):
        logger.critical(f'===== {len(setting.target_proteins) - len(rad_list)}개 단백질 RADetect 생성 실패, 프로그램 종료 =====')
        return

    rad : RADetect = rad_list[0]

    logger.info('===RESULT========================')
    for target_rad in rad_list:
        logger.info('===== 유사도 연산 진입 =====')
        similarity = rad.similar(target_rad, setting.gamma)
        logger.info(f'{rad.prot_id}, {target_rad.prot_id} 유사도 : {similarity}')
    logger.info('===RESULT========================')

#endregion

if __name__ == '__main__':
    main()