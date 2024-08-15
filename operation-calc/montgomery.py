from condition import Condition as C
from mytype import Op
import fp2
import sympy

# See https://sike.org/files/SIDH-spec.pdf for reference


# Calculate x(2P),z(2P) from x(P),z(P),A24plus, C24
def xDBL_proj(c: C) -> Op:
    return (
        # 1 to 3
        2 * fp2.add()
        + fp2.square(c)
        # 4 to 6
        + fp2.square(c)
        + 2 * fp2.mul(c)
        # 7 to 9
        + 2 * fp2.add()
        + fp2.mul(c)
        # 10
        + fp2.mul(c)
    )


# Calculate x(2^eP), z(2^eP) from x(P), z(P), A24plus, C24
def xDBLe_proj(c: C, e: Op) -> Op:
    return e * xDBL_proj(c)


def xTPL_proj(c: C) -> Op:
    return (
        # 1 to 6
        4 * fp2.add()
        + 2 * fp2.square(c)
        # 7 to 12
        + 2 * fp2.add()
        + fp2.square(c)
        + 3 * fp2.mul(c)
        # 13 o 18
        + 3 * fp2.add()
        + fp2.square(c)
        + 2 * fp2.mul(c)
        # 19 to 22
        + fp2.add()
        + fp2.square(c)
        + 2 * fp2.mul(c)
    )


def xTPLe_proj(c: C, e: Op) -> Op:
    return e * xTPL_proj(c)


# Calculate x(P+Q),z(P+Q),x(2P),z(2P)
# from x(P),z(P),x(Q),z(Q),x(Q-P),z(Q-P),a24plus
def xDBLADD_proj(c: C) -> Op:
    return (
        # 1 to 7
        4 * fp2.add()
        + 2 * fp2.square(c)
        + fp2.mul(c)
        # 8 to 14
        + 3 * fp2.mul(c)
        + 4 * fp2.add()
        # 15 to 19
        + 3 * fp2.mul(c)
        + 2 * fp2.square(c)
    )


# Calculate x(P+Q),y(P+Q) from x(P),y(P),x(Q),y(Q),a,b
def xADD_affi(c: C) -> Op:
    return (
        # 9 to 13
        2 * fp2.add()
        + fp2.inv(c)
        + fp2.mul(c)
        + fp2.square(c)
        # 14 to 18
        + 3 * fp2.add()
        + 2 * fp2.mul(c)
        # 19 to 23
        + 2 * fp2.mul(c)
        + 3 * fp2.add()
        # 24 to 26
        + 2 * fp2.add()
    )


# Calculate x(P + mQ),z(P + mQ) from x(P),x(Q),x(P-Q),a24plus
def three_point_ladder(c: C, bit_length_of_m: Op) -> Op:
    return bit_length_of_m * xDBLADD_proj(c)


# Calculate x(mP) from x(P), a24plus
def point_ladder(c: C, bit_length_of_m: Op) -> Op:
    return three_point_ladder(c, bit_length_of_m=bit_length_of_m)


# Calculate j-invariants from (A:C)
def j_inv(c: C) -> Op:
    return (
        # 1 to 5
        fp2.square(c)
        + 3 * fp2.add()
        # 6 to 10
        + fp2.square(c)
        + fp2.mul(c)
        + 3 * fp2.add()
        # 11 to 15
        + fp2.square(c)
        + fp2.mul(c)
        + fp2.inv(c)
        + 2 * fp2.add()
        # 16
        + fp2.mul(c)
    )


# Calculate (A24plus, C24) of E/<P> from (A24plus, C24) of E
# and x(P),z(P) where P is an order 2 point on E.
def two_iso_curve(c: C, bit_length_of_p: Op) -> Op:
    prob_x_0 = 1 / 3
    return prob_x_0 * (
        # t0, K
        fp2.sqrt(c, bit_length_of_p=bit_length_of_p)
        + fp2.mul(c)
        + 2 * fp2.add()
        # A24plus, C24
        + 5 * fp2.add()
    ) + (1 - prob_x_0) * (2 * fp2.square(c) + fp2.add())


# Calculate x(phi(P)),z(phi(P)) of E/<P> from (A24plus, C24) of E
# and x(Q), z(Q) on E and K.
def two_iso_eval(c: C) -> Op:
    prob_x_0 = 1 / 3
    return prob_x_0 * (
        # t0, t1,t2
        2 * fp2.square(c)
        + 3 * fp2.add()
        # A24plus, C24
        + 3 * fp2.mul(c)
        + fp2.add()
    ) + (1 - prob_x_0) * (4 * fp2.mul(c) + 6 * fp2.add())


# Calculate (A24plus,C24) of E/<P4>, K1,K2,K3 from (A24plus,C24) of E and P4 of E of order 4
def four_iso_curve(c: C) -> Op:
    prob_x_one = 2 / 12
    prob_x_minus_one = 2 / 12

    return (
        prob_x_one * fp2.add()
        + prob_x_minus_one * 0
        + (1 - prob_x_one - prob_x_minus_one)
        * (
            # 1 to 3
            2 * fp2.add()
            + fp2.square(c)
            # 4 to 6
            + 2 * fp2.add()
            + fp2.square(c)
            # 7 to 9
            + fp2.add()
            + 2 * fp2.square(c)
        )
    )


# Calculate (x(phi(Q)),z(phi(Q))) from x(Q),z(Q) and (K1,K2,K3) of four_iso_curve
def four_iso_eval(c: C) -> Op:
    prob_x_one = 2 / 12
    prob_x_minus_one = 2 / 12

    return (
        prob_x_one
        * (
            # t0,t1,t2
            3 * fp2.add()
            + fp2.square(c)
            # X, Z
            + 5 * fp2.mul(c)
            + 2 * fp2.add()
        )
        + prob_x_minus_one
        * (
            # t0,t1,t2
            3 * fp2.add()
            + fp2.square(c)
            + fp2.mul(c)
            # X, Z
            + 3 * fp2.mul(c)
            + fp2.add()
        )
        + (1 - prob_x_one - prob_x_minus_one)
        * (
            # 1 to 4
            2 * fp2.add()
            + 2 * fp2.mul(c)
            # 5 to 8
            + 2 * fp2.add()
            + 2 * fp2.mul(c)
            # 9 to 12
            + 2 * fp2.square(c)
            + 2 * fp2.add()
            # 13, 14
            + 2 * fp2.mul(c)
        )
    )


# Calculate (A24plus,C24) of E/<R> and (x(phi(Q)),z(phi(Q))) where phi: E -> E/<R>
# from (A24plus, C24) of E and R, Q where the order of R is 2^k
def two_e_iso(c: C, k: Op, number_eval_point=1) -> Op:
    return (
        # multiplication
        k * (k - 2) / 4 * xDBL_proj(c)
        # isogeny
        + (k / 2 - 1) * four_iso_curve(c)
        # point evaluation for R
        + (k / 2 - 2) * four_iso_eval(c)
        # point evaluation for Q
        + (k / 2 - 1) * four_iso_eval(c) * number_eval_point
    )


# Calculate 2^k-isogeny by optimal strategy
def optimal_strategy(c: C, k: Op, number_eval_point=1) -> Op:
    n = k / 2
    return (
        (
            # Balanced Strategy
            1
            / (2 * sympy.log(2))
            * n
            * sympy.log(n)
            # Mean Calculation cost of [4]P and 4-isogeny evaluation
            * (four_iso_eval(c) + xDBLe_proj(c, 2))
            / 2
        )
        # Bottom isogeny computation
        + (n - 1) * (four_iso_curve(c) + four_iso_eval(c))
        # Evaluate optional point
        + (n - 1) * four_iso_eval(c) * number_eval_point
    )


def three_iso_curve(c: C) -> Op:
    return (
        # 1 to 5
        3 * fp2.add()
        + 2 * fp2.square(c)
        # 6 to 10
        + 4 * fp2.add()
        + fp2.square(c)
        # 11 to 15
        + 4 * fp2.add()
        + fp2.mul(c)
        # 16 to 19
        + 2 * fp2.add()
        + fp2.mul(c)
    )


def three_iso_eval(c: C) -> Op:
    return (
        # 1 to 3
        2 * fp2.add()
        + fp2.mul(c)
        # 4 to 6
        + 2 * fp2.add()
        + fp2.mul(c)
        # 7 to 9
        + 2 * fp2.square(c)
        + fp2.mul(c)
        # 10
        + fp2.mul(c)
    )


def three_e_iso(c: C, e: Op, number_eval_point=1) -> Op:
    return (
        e * (e - 1) / 2 * xTPL_proj(c)
        + e * three_iso_curve(c)
        + (e - 1) * three_iso_eval(c)
        + e * three_iso_eval(c) * number_eval_point
    )