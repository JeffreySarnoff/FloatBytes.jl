```@meta
CurrentModule = FloatBytes
DocTestSetup = :(using FloatBytes)
```

# Projection: rounding and saturation

Every numeric operation follows one pipeline: decode the operands,
compute exactly over the closed extended reals, then apply **one final
projection** into the result format. A projection specification ρ is a
pair of a rounding mode and a saturation mode:

```julia
Projection(NearestTiesToEven(), SatNone())   # the package default
```

Internally projection is *round to precision*, then *saturate*, then
*encode* — never merged. Rounding happens first with unbounded exponent
range, so values between the maximum finite value and half an ulp above
it round back into range (draft §4.7.3 NOTE 1):

```jldoctest proj
julia> using FloatBytes

julia> quantize(Binary8p4se, 231.0)          # Mhi = 224, half-ulp above is 240
Binary8p4se(0x7e ↦ 224.0)

julia> quantize(Binary8p4se, 240.0)          # tie: rounds up, overflows → +Inf
Binary8p4se(0x7f ↦ Inf)
```

## Rounding modes

| mode | rule |
|:--|:--|
| [`NearestTiesToEven`](@ref) | nearest; ties to the even code (`P == 1`: evenness of the code point) |
| [`NearestTiesToAway`](@ref) | nearest; ties away from zero |
| [`TowardPositive`](@ref) | directed toward `+∞` |
| [`TowardNegative`](@ref) | directed toward `−∞` |
| [`TowardZero`](@ref) | truncation |
| [`ToOdd`](@ref) | truncate, force the kept code odd when inexact |
| [`StochasticA`](@ref)`(N)` | away iff `⌊ν·2^N⌋ + R ≥ 2^N` |
| [`StochasticB`](@ref)`(N)` | away iff `⌊ν·2^(N+1)⌋ + 2R + 1 ≥ 2^(N+1)` |
| [`StochasticC`](@ref)`(N)` | away iff `RNITE(ν·2^N) + R ≥ 2^N` |

`ν` is the exact fraction discarded below the ulp and `R < 2^N` the
supplied entropy value. `N ∈ 1:63` is validated at mode construction;
this build implements `N ≤ 32` (a larger `N` throws before any
computation — see [Stochastic rounding](stochastic.md)).

## Saturation modes

| mode | out-of-range behavior |
|:--|:--|
| [`SatFinite`](@ref) | clamp everything, including infinities, to the extremal finite values |
| [`SatPropagate`](@ref) | representable infinities propagate; finite overflow clamps |
| [`SatNone`](@ref) | draft §4.7.5 default: `±Inf` where representable, NaN where not (finite domains; negative results in unsigned formats) — but a rounding direction pointing inward forces the extremal finite value |

```jldoctest proj
julia> quantize(Binary8p4sf, 1.0e6)                       # finite domain, SatNone
Binary8p4sf(0x80 ↦ NaN)

julia> quantize(Binary8p4sf, 1.0e6, Projection(NearestTiesToEven(), SatFinite()))
Binary8p4sf(0x7f ↦ 240.0)

julia> quantize(Binary8p4se, 1.0e6, Projection(TowardZero(), SatNone()))
Binary8p4se(0x7e ↦ 224.0)
```

## Scoped override for Base operators

The Base operator surface (`+`, `*`, `sqrt`, …) reads a scoped default;
[`with_projection`](@ref) overrides it for a dynamic extent, with no
mutable global state:

```jldoctest proj
julia> with_projection(() -> Binary8p4se(1.0) / Binary8p4se(3.0),
                       Projection(TowardPositive(), SatNone()))
Binary8p4se(0x33 ↦ 0.34375)

julia> with_projection(() -> Binary8p4se(1.0) / Binary8p4se(3.0),
                       Projection(TowardNegative(), SatNone()))
Binary8p4se(0x32 ↦ 0.3125)
```

The explicit [`apply`](@ref)`(op, Fr, ρ, xs...)` form is the primary
interface; the scope only affects the Base-operator convenience surface.
