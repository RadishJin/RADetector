# RADetector

---
RADetector is a lightweight pipeline designed for benchmarking protein structural evaluation metrics:
RMSD, lDDT, TM-score

Made by Mu Jin Kim (BS student @ Kyungpook National University, Department of Biotechnology)
Contact : radishj24@gmail.com
Managed under GitFlow branching model
---


 - Versions
    Python 3.11
    ...


 - Directories
    RADetector/

        .vscode/
            settings.json

        raw_dataset/        # Raw dataset from PDB [Ignored by Git]
        test_dataset/       # Noise-injected decoy datasets [Ignored by Git]

        result/             # Final evaluation summary
        src/
            gathering.sh    # Bash script for fetching datasets
            parsing.py      # Data parsing and tensor conversion
            preparing.py    # Noise injection and decoy generation
            testing.py      # Metric calculation (TM-score, lDDT, RMSD)

        .gitignore
        LICENSE
        README.md


 - Brief Explanation
    Algorithm
    1. 데이터 전처리 (정답 데이터셋 만들기)
        1. Bash, Protein Dataset 가져오기
        2. Python, Data Parsing
        3. Python, PyTorch Tensor로 저장

    2. 테스트 데이터셋 만들기
        1. 저장한 PyTorch Tensor 데이터 가져오기, 외부 실험 데이터 가져오기
        2. 정답 데이터에 특정 강도 노이즈 주입
        3. 테스트 데이터로 저장

    3. 결과 정리하기
        1. 정답 데이터, 테스트 데이터 가져오기
        2. TM-score, lDDT, RMSD 산출
        3. 결과 요약 출력

    Target Protein
	    Standard Rigid Target - (1CRN - Crambin)
	    Domain-Flexible Target - (1CLL - Calmodulin)
	    Antibody - (5DK3 - Pembrolizumab Fab)


