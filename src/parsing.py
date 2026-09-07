import gzip
import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc

# gzip으로 압축파일 간단하게 열고, biotite 라이브러리 이용해서 객체로 바로 읽어오기
with gzip.open("raw_dataset/1CLL.cif.gz", "rt", encoding = "utf-8") as f:
    raw = pdbx.CIFFile.read(f)
# print(raw)


# pdbx.get_structrue로 텍스트데이터를 AtomArray로 파싱
atoms = pdbx.get_structure(raw, model= 1)
# print(atoms)

# Biotite.structure 이용 residue만 남기는 불리언 마스크 생성
residues = struc.filter_amino_acids(atoms)
# print(residues)

# 불리언 마스크 통해 유효 residue의 원자들만 남기기
chain = atoms[residues]
# print(chain)

# Canonicalization



# 데이터 중 원자 데이터만 가져오기


# .cif 자료구조 (일부)

# loop_
# _atom_site.group_PDB      # 1번 컬럼: ATOM 여부
# _atom_site.id             # 2번 컬럼: 원소 번호
# _atom_site.type_symbol    # 3번 컬럼: 원소 기호 (C, N, O 등)
# _atom_site.label_atom_id  # 4번 컬럼: 원자 이름 (CA, CB 등)
# _atom_site.label_comp_id  # 5번 컬럼: 잔기 이름 (ALA, GLY 등)
# _atom_site.label_seq_id   # 6번 컬럼: 잔기 번호
# _atom_site.Cartn_x        # 7번 컬럼: X 좌표
# _atom_site.Cartn_y        # 8번 컬럼: Y 좌표
# _atom_site.Cartn_z        # 9번 컬럼: Z 좌표

