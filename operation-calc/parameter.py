from sympy import symbols, simplify
from mytype import Op

## The cost for calculating multiplication, square, add, inversion in Fp
M, S, A, I = symbols("M S A I")


def M_part(op: Op) -> Op:
    return simplify(op.subs([(S, 0), (A, 0), (I, 0)]) / M)


def S_part(op: Op) -> Op:
    return simplify(op.subs([(M, 0), (A, 0), (I, 0)]) / S)


def A_part(op: Op) -> Op:
    return simplify(op.subs([(M, 0), (S, 0), (I, 0)]) / A)


def I_part(op: Op) -> Op:
    return simplify(op.subs([(M, 0), (S, 0), (A, 0)]) / I)


def to_str(op: Op) -> str:
    msai_repr = ""
    m_part = M_part(op)
    if m_part != 0:
        msai_repr += f"({m_part.evalf(6)})M"

    s_part = S_part(op)
    if s_part != 0:
        if msai_repr != "":
            msai_repr += " + "
        msai_repr += f"({s_part.evalf(6)})S"

    a_part = A_part(op)
    if a_part != 0:
        if msai_repr != "":
            msai_repr += " + "
        msai_repr += f"({a_part.evalf(6)})A"

    i_part = I_part(op)
    if i_part != 0:
        if msai_repr != "":
            msai_repr += " + "
        msai_repr += f"({i_part.evalf(6)})I"

    return msai_repr


def print_MSAI(op: Op):
    print(f"M: {M_part(op)}")
    print(f"S: {S_part(op)}")
    print(f"A: {A_part(op)}")
    print(f"I: {I_part(op)}")