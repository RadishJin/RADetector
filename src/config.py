from __future__ import annotations

# Custom Mother Exception
# =========================================
class PipelineError(Exception):
    """ 커스텀 부모 클래스 에러 """
# =========================================


# Custom Daughter Exceptions
# =========================================
class NotFoundError(PipelineError):
    """ 데이터가 존재하지 않을 때 발생하는 예외 """

class NotValidError(PipelineError):
    """ 데이터가 형태에 맞지 않을 때 발생하는 예외 """

class EmbedFailError(PipelineError):
    """ 분자 3D 임베딩 실패 """
# =========================================

    