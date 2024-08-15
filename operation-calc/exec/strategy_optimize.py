from strategy_optimization import StrategyOptimizer
import montgomery as mont
from condition import Condition
from mytype import Op
import fp


def n_to_size_p(n: int) -> int:
    if n <= 256:
        return 256
    elif n <= 512:
        return 512
    elif n <= 1024:
        return 1024
    elif n <= 1536:
        return 1536
    else:
        assert False


def main():
    c = Condition(is_p_equal_3_mod_4=True)
    mul_by_2 = mont.xDBL_proj(c)
    mul_by_4 = mul_by_2 * 2
    point_eval_4isog = mont.four_iso_eval(c)
    curve_eval_4isog = mont.four_iso_curve(c)

    reduce_maps: dict[int, dict[str, float]] = {
        256: {"A": 0.24, "S": 1.02, "I": 105.39},
        512: {"A": 0.13, "S": 1.00, "I": 118.69},
        1024: {"A": 0.050, "S": 1.00, "I": 101.73},
        1536: {"A": 0.029, "S": 1.00, "I": 86.93},
    }

    print("n: ")
    n = int(input())
    reduce_map = reduce_maps[n_to_size_p(n)]

    mul_by_2 = fp.reduce(mul_by_2, reduce_map) / fp.mul()
    mul_by_4 = fp.reduce(mul_by_4, reduce_map) / fp.mul()
    point_eval_4isog = fp.reduce(point_eval_4isog, reduce_map) / fp.mul()
    curve_eval_4isog = fp.reduce(curve_eval_4isog, reduce_map) / fp.mul()

    print("Cost to compute some functions")
    print("mul_by_2: ", mul_by_2)
    print("mul_by_4: ", mul_by_4)
    print("point_eval_4isog: ", point_eval_4isog)
    print("curve_eval_4isog: ", curve_eval_4isog)

    optimizer = StrategyOptimizer(
        mul_by_4, mul_by_2, point_eval_4isog, curve_eval_4isog
    )

    cost, strategy = optimizer.optimize_backtrack(n)
    print(f"Cost: {cost}M")
    print(f"Strategy: {strategy}")


if __name__ == "__main__":
    main()
