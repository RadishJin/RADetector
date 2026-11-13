from __future__ import annotations

class PipelineError(Exception):
    """ 파이프라인 에러 전체의 모예외 """

class NotFoundError(PipelineError):
    """ 데이터가 존재하지 않을 때 발생하는 예외 """

class EmbedFailError(PipelineError):
    """ 분자 임베딩 실패 시 발생하는 예외 """

class NotValidError(PipelineError):
    """ 데이터가 형태에 맞지 않을 때 발생하는 예외 """
