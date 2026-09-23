# RADetector
Identifying functional similarity between proteins through their molecular-interaction patterns.
Author : MU JIN KIM, BS student @ kyungpook national univ.
Contact : radishj24@gmail.con


Main Protein : EGFR (Epidermal Growth Factor Receptor) CHEMBL203

Comparison A : HER2 (Human Epidermal Growth Factor Receptor 2) CHEMBL1824
 > Same Kinase Family, Similar RADetect
Comparison B : CSF1R (Macrophage colony-stimulating factor 1 receptor) CHEMBL1844
 > Same Kinase Family, Different RADetect
Comparison C : CDK2 (Cyclin-dependent kinase 2) CHEMBL301
 > Different Family, Different RADetect
Comparison D : CA2 (Carbonic Anhydrase 2) CHEMBL205
 > Different Family, Simular RADetect


Overview

1. target protein의 molecule interaction data parsing (bash, python)
    i. API로 데이터셋 긁어오기
    ii. 결합 분자 리스트 만들기

2. target protein의 molecular interaction 데이터로 일종의 지문 도출 (python)
    ii. 가로축은 Van der Waals volume, 세로축은 pKa로 ndarray 배열 생성
    iii. Gaussian Kernal Density Estimation 이용 2D 이미지로 변환
    v. target protein에 할당

3. ndarray 유사도 평가 (python)
    i. gaussian kde를 cosine simularity 비교
    ii. 문서화

4. 실제 알려진 기능 유사도와 RADetect 결과물 비교 (결과 작성)


