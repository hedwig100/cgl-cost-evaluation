# CGL cost evaluation

This repository contains codes on "A faster variant of CGL hash function via efficient backtracking checks" at Information Security Conference 2024.

This repository consists of two packages. First package is made of C++ and protocols are implemented. Second package is made of Python and it is to calculate the number of operations or optimize isogeny computation strategy.

## Cost evaluation of CGL hash functions

After you install the following requirements, run `./do run count-cgl`. You will get the result of performance evaluation.

NOTE: DPB-CGL is **NOT** collision-resistant, it is used just for comparison.

```
```

## Collision finding when the backtracking check is lost

We show the collision can be found when the backtracking check is lost as in Appendix in our paper.
Run `./do run find_collision`. This command sometimes fails. If it succeeds, you will get a collision pair.

```
```

## Requirements [C++]

We need `libboost-dev` as a multiprecision computation library.

```
sudo apt-get install libboost-dev
```

## Execution [C++]

First, build

```
./do build
```

Then run test

```
./do test
```

Finally, run any `src/<exec>.cpp`

```
./do run <exec>
```

## Requirements [Python]

Install poetry.

## Exec [Python]

First, install libraries and run a sample python code.

```
poetry install
```

Then, run scirpts.

```
./do python strategy_optimize
```