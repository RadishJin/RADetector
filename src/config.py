from __future__ import annotations


# =========================================
# Custom Exceptions
class NotFoundError(Exception):
    """ 데이터가 존재하지 않을 때 발생하는 예외 """

class NotValidError(Exception):
    """ 데이터가 형태에 맞지 않을 때 발생하는 예외 """

class EmbedFailError(Exception):
    """ 분자 3D 임베딩 실패 """
# =========================================

    