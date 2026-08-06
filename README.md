# FloatBytes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JeffreySarnoff.github.io/FloatBytes.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JeffreySarnoff.github.io/FloatBytes.jl/dev/)
[![Build Status](https://github.com/JeffreySarnoff/FloatBytes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JeffreySarnoff/FloatBytes.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JeffreySarnoff/FloatBytes.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JeffreySarnoff/FloatBytes.jl)

IEEE P3109 (working draft D1) byte-sized floating point for Julia 1.12+:
all 120 `Binary{K}p{P}{s|u}{e|f}` formats for `K ∈ 3:8`, with exact
integer-window arithmetic, draft-literal special-value semantics, every
draft rounding mode (including stochastic), and saturation as a
first-class projection component.

This is an implementation of a **draft** standard; it does not claim
P3109 conformance. Design: `docs/design/ImplementationPlan.md`; current
coverage: `docs/design/ImplementationStatus.md`.

```julia
julia> using FloatBytes

julia> x = Binary8p4se(1.5)          # projecting constructor (NTE, SatNone)
Binary8p4se(0x44 ↦ 1.5)

julia> x + Binary8p4se(2.25)         # Base operators, same format
Binary8p4se(0x4f ↦ 3.75)

julia> apply(Divide, Binary8p4se, Projection(TowardZero(), SatFinite()), x, x + x)
Binary8p4se(0x38 ↦ 0.5)

julia> quantize(Binary5p2ue, 1.06,   # stochastic rounding: entropy is input
                Projection(StochasticA(8), SatNone());
                entropy = IndexedEntropy(42), index = 0)
Binary5p2ue(0x10 ↦ 1.0)

julia> one(Binary8p4se) / zero(Binary8p4se)   # draft semantics: x/0 = NaN
Binary8p4se(0x80 ↦ NaN)

julia> sort([Binary3p1se(c) for c in rand(Float64, 4)])  # NaN-first total order
```

Key departures from IEEE-754 habits (all normative in the draft): one
NaN and one zero per format, `x/0 = NaN`, `Recip(±Inf) = 0`,
`NextGreaterThan(+Inf) = NaN`, NaN sorts *first* in the total order, and
finite-domain formats saturate per the projection's saturation mode.
