```@meta
CurrentModule = FloatBytes
DocTestSetup = :(using FloatBytes)
```

# Operations

The public seam is [`apply`](@ref):

```julia
apply(op, Fr, ρ, xs...; entropy = nothing, index = 0) -> Fr
```

`op` is an operation singleton ([`Add`](@ref), [`FMA`](@ref),
[`Sqrt`](@ref), …), `Fr` the result format, `ρ` the projection.
Operands are `FloatByte` values and may have **different formats** —
result format and operand formats are all explicit. Special values
follow the draft's pattern tables in their normative order; finite
arguments are evaluated **exactly** (integer window arithmetic) and
projected once.

```jldoctest ops
julia> using FloatBytes

julia> ρ = DEFAULT_PROJECTION;

julia> apply(FMA, Binary8p4se, ρ, Binary8p4se(1.5), Binary8p4se(2.25), Binary8p4se(2.25))
Binary8p4se(0x53 ↦ 5.5)

julia> apply(Add, Binary8p3se, ρ, Binary8p4se(1.5), Binary5p2ue(0.5))  # mixed formats
Binary8p3se(0x44 ↦ 2.0)
```

## Operation families

- **Arithmetic** (draft §4.10): [`Add`](@ref), [`Subtract`](@ref),
  [`Multiply`](@ref), [`Divide`](@ref), [`FMA`](@ref), [`FAA`](@ref),
  [`Sqrt`](@ref), [`Recip`](@ref), [`RSqrt`](@ref), [`Hypot`](@ref),
  [`Abs`](@ref), [`Negate`](@ref), [`CopySign`](@ref),
  [`Convert`](@ref).
- **Extrema** (§4.11): [`Minimum`](@ref)/[`Maximum`](@ref) and the
  `Number`, `Magnitude`, `MagnitudeNumber`, and `Finite` variants, plus
  [`Clamp`](@ref). Ten min/max variants total, each with its own NaN and
  infinity rows.
- **Comparisons and predicates** (§4.12–4.13): [`compare_less`](@ref)
  and friends, [`total_order`](@ref), `isnan`/`isinf`/`isfinite`/
  `issubnormal`/`iszero`/`isone`/`signbit`, and the classifier
  [`floatclass`](@ref). These return Booleans or an enum — no
  projection, no entropy.
- **Neighbours** (§4.16): [`nextgreaterthan`](@ref),
  [`nextlessthan`](@ref) — format-preserving code-point walks with
  explicit terminal rows.
- **Transcendentals** (§4.10.9–.16): `Exp` … `ArcTan2Pi` — see below.

## Fused operations project once

[`FMA`](@ref) and [`FAA`](@ref) assemble the complete expression exactly
and round a single time. FAA is exactly associative — bit-identical for
every operand order:

```jldoctest ops
julia> x = Binary8p4se(224.0); y = Binary8p4se(-224.0); z = Binary8p4se(0.001953125);

julia> apply(FAA, Binary8p4se, ρ, x, y, z) == apply(FAA, Binary8p4se, ρ, z, x, y)
true

julia> Float64(apply(FAA, Binary8p4se, ρ, x, y, z))   # exact cancellation survives
0.001953125
```

## Transcendentals: correct rounding, or an error

Each of the 25 unary transcendentals plus `ArcTan2`/`ArcTan2Pi` runs a
rigorous enclosure loop: directed MPFR endpoints at doubling precision
until the projection is decided. There is **no undeclared
approximation** — if the loop cannot decide at its precision ceiling it
throws [`InconclusiveProjectionError`](@ref) rather than writing a
guess. Exact points (`Exp(0) = 1`, `Log(1) = 0`, `SinPi` at
half-integers, `TanPi(±1/4) = ±1`, `Tanh(±Inf) = ±1`, …) are exact, and
asymptotes keep their sidedness (`Tanh` of a large finite value is
`1 − ε`, never `1`, so directed rounding toward zero stays below one):

```jldoctest ops
julia> Float64(apply(TanPi, Binary8p4se, ρ, Binary8p4se(0.25)))
1.0

julia> apply(Tanh, Binary8p4se, Projection(TowardZero(), SatNone()), maxfinite(Binary8p4se))
Binary8p4se(0x3f ↦ 0.9375)
```
