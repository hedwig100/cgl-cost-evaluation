import parameter as param
from condition import Condition as C
from mytype import Op


def reduce(
    op: Op, reduce_map: dict[str, float] = {"S": 0.8, "I": 10, "A": 0.004}
) -> Op:
    return op.subs(
        [
            (param.S, reduce_map["S"] * param.M),
            (param.I, reduce_map["I"] * param.M),
            (param.A, reduce_map["A"] * param.M),
        ]
    )


def mul():
    return param.M


def add():
    return param.A


def square():
    return param.S


def inv():
    return param.I


def div():
    return inv() + mul()


# x^n in Fp
def exp(c: C, bit_length_of_n: Op, stand_bit_of_n: Op) -> Op:
    return bit_length_of_n * param.S + stand_bit_of_n * param.M


# Check if x is the QNR in Fp
def chi(c: C, bit_length_of_p: Op) -> Op:
    return exp(
        c, bit_length_of_n=bit_length_of_p - 1, stand_bit_of_n=(bit_length_of_p - 1) / 2
    )


# Calculate √x of x in Fp
def sqrt(c: C, bit_length_of_p: Op):
    if c.is_p_equal_3_mod_4:
        return exp(
            c,
            bit_length_of_n=bit_length_of_p - 2,
            stand_bit_of_n=(bit_length_of_p - 2) / 2,
        )
    # TODO: add tonelli-shanks operation count
    raise NotImplementedError