# cocks-pinch

Construct pairing-friendly elliptic curves using the Cocks-Pinch method with Complex Multiplication (CM).

Given a prime subgroup order `r` and target embedding degree `k`, finds a prime `p` and trace `t` such that `r | p + 1 - t` and `p = (t^2 + D*y^2) / 4` for a class-number-1 CM discriminant `D`.

## Build

Requires Zig >= 0.15.2.

```
zig build
```

## Run

```
zig build run
```

Generates `curve.txt` and `curve.json` with full curve parameters.

Default configuration: secp256k1 subgroup order, embedding degree k=6, targeting ~128-bit security (F_p^6 ~ 2^3156).

## Test

```
zig build test
```

## Library

Core functions in `src/root.zig`:

- `cocksPinchCM` -- main entry point, searches over CM discriminants (D = 3, 4, 7, 8, 11, 19, 43, 67, 163)
- `cocksPinch` -- basic Cocks-Pinch without CM restriction
- `constructCurveFromJ` -- derives curve coefficients `a`, `b` from j-invariant
- `powMod`, `sqrtMod`, `legendreSymbol`, `modInverse` -- modular arithmetic primitives
- `isPrime` -- Miller-Rabin primality test
- `findKthRootOfUnity` -- primitive k-th root of unity mod r

Import as a Zig module:

```zig
const cp = @import("cocks_pinch");
```
