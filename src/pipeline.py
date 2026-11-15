#region MODULE & LIBRARY
from __future__ import annotations

# Standard Libraries
import logging

# Custom Modules & Classes
from io_parser import io_prot_bindmol, parse_mol_query
from domain import get_known_props, Molecule, RADetect

from config import Config, Setting
from exceptions import DataError, PipelineError

#endregion
#region GLOBAL VARIABLE

logger = logging.getLogger(__name__)

#endregion
#region MAIN

def radetector(prot_id : str, setting : Setting, config : Config) -> RADetect:
    logger.info(f'===== Protein : {prot_id} - RADetect 생성 진입 =====')

    try:
        mol_raw = io_prot_bindmol(prot_id)
        logger.debug('Binding Assay 리스트 파싱 성공')
    except DataError as e:
        logger.error(f'Protein : {prot_id} - Binding Assay 수신 실패: {e}')
        raise PipelineError(e)

    mol_id_list = list(set([raw['molecule_chembl_id'] for raw in mol_raw]))
    logger.info(f'Protein : {prot_id} - 총 {len(mol_raw)} 중 유니크 분자 {len(mol_id_list)}개 추출, 물성치 조회 시작')
    
    try:
        known_props = get_known_props(mol_id_list, setting.batch_size, setting.max_retries)
        logger.debug('알려진 Molecular Property 리스트 수신 성공')
    except DataError as e:
        logger.error(f'Protein : {prot_id} - Molecular Property 수신 실패: {e}')
        raise PipelineError(e)

    seen_ids = set()
    for raw in mol_raw:
        data : tuple = parse_mol_query(raw)

        if (data[0] is None) or (data[1] is None) or (data[2] is None):
            continue

        mol_id = data[0]

        if mol_id in seen_ids:
            continue
        seen_ids.add(mol_id)

        try:
            mol = Molecule.extended(mol_id, data[1], data[2], setting.maxseed, setting.maxiter)
        except DataError as e:
            config.omits[prot_id].append(mol_id)
            logger.error(f'실패한 [Protein : {prot_id}, Molecule : {mol_id}] - {e}')
            continue
        
        props = known_props.get(mol.id)
        if props:
            mol.properties.update(props)
            config.res_mols.append(mol)
        else:
            logger.warning(f'물성치 누락 [Protein : {prot_id}, Molecule : {mol.id}] - 알려진 Property 없음')

    logger.info(f'Protein : {prot_id} - {len(config.res_mols)}/{len(mol_id_list)}개 분자 파싱 성공')

    config.res_prots[prot_id] = list(config.res_mols)
    config.res_mols.clear()

    return RADetect.extended(prot_id, config.res_prots[prot_id])

#endregion