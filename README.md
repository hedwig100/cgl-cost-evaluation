# CGL cost evaluation

This repository contains codes on "A faster variant of CGL hash function via efficient backtracking checks" at Information Security Conference 2024.

This repository consists of two packages. First package is made of C++ and protocols are implemented. Second package is made of Python and it is to calculate the number of operations or optimize isogeny computation strategy.

## Cost evaluation of CGL hash functions

After you install the following requirements, run `./do run count-cgl`. You will get the result of performance evaluation.

NOTE: DPB-CGL is **NOT** collision-resistant, it is used just for comparison.

```
...
Consolidate compiler generated dependencies of target count-cgl
[100%] Built target count-cgl
+ ./build/src/count-cgl
Size of p (bit):256
n:240
P: 94482146789024102024466833852227567633874701385423723785552891614122842521599
e: 240
  name: ModularPolynomial-CGL
  Count:   M_COUNTER: 279.57
  S_COUNTER: 758.615
  A_COUNTER: 24.445
  I_COUNTER: 0.985
  name: Yoshida-Takashima-CGL
  Count:   M_COUNTER: 273.095
  S_COUNTER: 762.45
  A_COUNTER: 12.095
  I_COUNTER: 0.995
  name: Hashimoto-Nuida-CGL
  Count:   M_COUNTER: 383.41
  S_COUNTER: 510.88
  A_COUNTER: 7.835
  I_COUNTER: 0.505
  name: Radical-CGL
  Count:   M_COUNTER: 388.73
  S_COUNTER: 510.925
  A_COUNTER: 9.165
  I_COUNTER: 0.505
  name: DPB-CGL
  Count:   M_COUNTER: 165.076
  S_COUNTER: 3.43542
  A_COUNTER: 176.323
  I_COUNTER: 0.004375
  name: DPB-CGL-prevent-backtrack
  Count:   M_COUNTER: 205.565
  S_COUNTER: 6.41
  A_COUNTER: 215.793
  I_COUNTER: 0.00833333
  name: My-CGL
  Count:   M_COUNTER: 199.82
  S_COUNTER: 5.24104
  A_COUNTER: 211.639
  I_COUNTER: 0.00666667
...
```

## Collision finding when the backtracking check is lost

We show the collision can be found when the backtracking check is lost as in Appendix in our paper.
Run `./do run find_collision`. This command sometimes fails. If it succeeds, you will get a collision pair.

```
...
Consolidate compiler generated dependencies of target find_collision
[100%] Built target find_collision
+ ./build/src/find_collision
Found collision!!
P:112200089154621740081528141189677554510592709903292430697596358600685233635327
n: 240
m0: 101001001100110110001110100011000001110100111100111101000111000001111011101000000111101100110110011101001000101101100101011011111111000111101010010011111100010010010001001111011000101111011011111101111011000101110101101100110001110011110010
m0': 001110010110110001000110100101011111101110011100110010000011011111101110001101000000110111101010111010011010000010010100011011111111000111101010010011111100010010010001001111011000101111011011111101111011000101110101101100110001110011110010
m1': 110111000011011001111001011010110110101011110010000111101011011010000101010000111110001111001000111000101110000000110100011101000101100110100011110001100111000110111001111111110000101001101111101000010101011111001100110001101011110001101111
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