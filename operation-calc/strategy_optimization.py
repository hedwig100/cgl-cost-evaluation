INF: float = 1e18


class StrategyOptimizer:
    def __init__(
        self,
        mul_by_4: float,
        mul_by_2: float,
        point_eval_4isog: float,
        curve_eval_4isog: float,
    ) -> None:
        self.mul_by_4 = mul_by_4
        self.mul_by_2 = mul_by_2
        self.point_eval_4isog = point_eval_4isog
        self.curve_eval_4isog = curve_eval_4isog

    def optimize_backtrack(self, n: int):
        if n % 2 == 0:
            return self.optimize_backtrack_even(n)
        else:
            return self.optimize_backtrack_odd(n)

    def optimize(self, n: int):
        if n % 2 == 0:
            return self.optimize_even(n)
        print("Odd n is not supported yet.")
        return [], []

    def _optimize_even(self, n: int):
        k = n // 2

        # dp0[i] = min cost to compute 4^i-isogeny
        dp0 = [INF] * (k + 1)
        from_where = [0] * (k + 1)

        dp0[1] = self.curve_eval_4isog
        for i in range(2, k + 1):
            for j in range(1, i):
                if (
                    dp0[j]
                    + dp0[i - j]
                    + (i - j) * self.mul_by_4
                    + j * self.point_eval_4isog
                    < dp0[i]
                ):
                    dp0[i] = (
                        dp0[j]
                        + dp0[i - j]
                        + (i - j) * self.mul_by_4
                        + j * self.point_eval_4isog
                    )
                    from_where[i] = j
        return dp0, from_where

    def _retrieve_strategy_even(self, from_where: list[int], k: int):
        if k == 1:
            return []
        strategy = [(k - from_where[k])]
        strategy.extend(self._retrieve_strategy_even(from_where, from_where[k]))
        strategy.extend(self._retrieve_strategy_even(from_where, k - from_where[k]))
        return strategy

    def optimize_even(self, n: int):
        dp0, from_where = self._optimize_even(n)
        k = n // 2
        return dp0[k], self._retrieve_strategy_even(from_where, k)

    def optimize_backtrack_even(self, n: int):
        dp0, from_where = self._optimize_even(n)
        k = n // 2
        return dp0[k], [2 * a for a in self._retrieve_strategy_even(from_where, k)]

    def _retrieve_strategy_odd(
        self, from_where0: list[int], from_where1: list[int], k: int
    ):
        if k == 0:
            return []
        strategy = [2 * (k - from_where1[k]) + 1]
        strategy.extend(
            [2 * a for a in self._retrieve_strategy_even(from_where0, from_where1[k])]
        )
        strategy.extend(
            self._retrieve_strategy_odd(from_where0, from_where1, k - from_where1[k])
        )
        return strategy

    def optimize_backtrack_odd(self, n: int):
        dp0, from_where0 = self._optimize_even(n - 1)
        k = (n - 1) // 2

        # dp1[i] = min cost to compute 2\cdot 4^i = 2^{2i + 1}-isogeny
        dp1 = [INF] * (k + 1)
        from_where1 = [0] * (k + 1)
        dp1[0] = 0
        for i in range(1, k + 1):
            for j in range(1, i + 1):
                if (
                    dp0[j]
                    + dp1[i - j]
                    + self.mul_by_2
                    + (i - j) * self.mul_by_4
                    + j * self.point_eval_4isog
                    < dp1[i]
                ):
                    dp1[i] = (
                        dp0[j]
                        + dp1[i - j]
                        + self.mul_by_2
                        + (i - j) * self.mul_by_4
                        + j * self.point_eval_4isog
                    )
                    from_where1[i] = j

        return dp1[k], self._retrieve_strategy_odd(from_where0, from_where1, k)
