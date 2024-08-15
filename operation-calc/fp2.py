import parameter as param
from condition import Condition as C
from mytype import Op
import fp


def mul(c: C) -> Op:
    if c.is_p_equal_3_mod_4:
        # (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        return 4 * param.M + 2 * param.A
    else:
        # (a + bi)(c + di) = (ac + βbd) + (ad + bc)i
        return 5 * param.M + 2 * param.A


def square(c: C) -> Op:
    if c.is_p_equal_3_mod_4:
        # (a + bi)^2 = (a + b)(a - b) + 2abi
        return 2 * param.M + 3 * param.A
    else:
        # (a + bi)^2 = (a^2 + βb^2) + 2abi
        return 2 * param.M + 2 * param.S + 2 * param.A


def add() -> Op:
    # (a + bi) + (c + di) = (a + c) + (b + d)i
    return 2 * param.A


def inv(c: C) -> Op:
    if c.is_p_equal_3_mod_4:
        # 1/(a + bi) = (a - bi)/(a^2 + b^2)
        return 2 * param.M + 2 * param.S + 2 * param.A + param.I
    else:
        # 1/(a + bi) = (a - bi)/(a^2 - βb^2)
        return 3 * param.M + 2 * param.S + 2 * param.A + param.I


def div(c: C) -> Op:
    return inv(c) + mul(c)


# x^n in Fp2
def exp(c: C, bit_length_of_n: Op, stand_bit_of_n: Op) -> Op:
    return bit_length_of_n * square(c) + stand_bit_of_n * mul(c)


# Calcaulte N(a + bi) = (a + bi)^((p+1)/2) = a^2 - beta b^2
# where i^2 = beta
def norm(c: C) -> Op:
    if c.is_p_equal_3_mod_4:
        return param.A + 2 * param.S
    else:
        return param.A + 2 * param.S + param.M


# Calculate x^{(p^2 - 1)/2} of x in Fp2
def chi(c: C, bit_length_of_p: Op) -> Op:
    return norm(c) + exp(
        c, bit_length_of_n=bit_length_of_p - 1, stand_bit_of_n=(bit_length_of_p - 1) / 2
    )


# Calculate √x of x in Fp2
def sqrt(c: C, bit_length_of_p: Op) -> Op:
    # TODO: implemente scott trick version
    return (
        norm(c)
        + fp.sqrt(c, bit_length_of_p)
        + fp.mul()
        + fp.add()
        + fp.chi(c, bit_length_of_p)
        + fp.add()
        + fp.sqrt(c, bit_length_of_p)
        + fp.mul()
        + fp.div()
    )