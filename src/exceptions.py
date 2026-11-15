from __future__ import annotations

class PipelineError(Exception):
    """ 파이프라인 최상위 예외 """

class DataError(PipelineError):
    """ 데이터 결손 시 발생하는 하위 예외 """

class NotFoundError(DataError):
    """ 데이터가 존재하지 않을 때 발생하는 최하위 예외 """

class EmbedFailError(DataError):
    """ 분자 임베딩 실패 시 발생하는 최하위 예외 """

class NotValidError(DataError):
    """ 데이터가 형태에 맞지 않을 때 발생하는 최하위 예외 """

class TimeOutError(DataError):
    """ API 요청 한도 초과 시 발생하는 최하위 예외 """