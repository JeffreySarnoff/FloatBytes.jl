```@meta
CurrentModule = FloatBytes
DocTestSetup = :(using FloatBytes)
```

# Formats

A format is a tuple `(K, P, Σ, Δ)`: bitwidth `K ∈ 3:8`, precision `P`,
signedness `Σ ∈ {Signed, Unsigned}`, domain `Δ ∈ {Finite, Extended}`.
The value type is the one-byte primitive
[`FloatByte`](@ref)`{K,P,S,D}`; users normally write the 120 exported
aliases, named per the draft:

> `Binary{K}p{P}{s|u}{e|f}` — e.g. `Binary8p4se` is 8 bits, precision 4,
> **s**igned, **e**xtended; `Binary5p2uf` is 5 bits, precision 2,
> **u**nsigned, **f**inite.

Signed formats allow `1 ≤ P ≤ K−1`; unsigned allow `1 ≤ P ≤ K`. Extended
formats have `+Inf` (and `−Inf` when signed) as code points; finite
formats have neither. Every format has exactly one NaN and one zero.

```jldoctest formats
julia> using FloatBytes

julia> bitwidth(Binary8p4se), precision(Binary8p4se), exponent_bias(Binary8p4se)
(8, 4, 8)

julia> signedness(Binary8p4se), domain(Binary8p4se)
(SignedFormat, ExtendedDomain)

julia> Float64.((maxfinite(Binary8p4se), minpositive(Binary8p4se)))
(224.0, 0.0009765625)
```

## Code identity vs. mathematical value

Two constructors, visibly different:

- [`fromcode`](@ref)`(F, code)` — **raw code identity**, range-checked.
  The only safe way to name a specific bit pattern.
- `F(value::Real)` / [`quantize`](@ref) — **mathematical construction**:
  projects a real value into the format (`F(value)` uses the package
  default projection; `quantize` takes an explicit ρ).

```jldoctest formats
julia> FloatBytes.fromcode(Binary8p4se, 0x44)
Binary8p4se(0x44 ↦ 1.5)

julia> FloatBytes.codepoint(Binary8p4se(1.5))
0x44

julia> quantize(Binary8p4se, 2/3, Projection(TowardNegative(), SatNone()))
Binary8p4se(0x3a ↦ 0.625)
```

`reinterpret` into a `FloatByte` is documented **unsafe**: it can
fabricate nonzero bits above `K`, violating the canonical-byte
invariant. `show` masks and flags such bytes instead of throwing.

## Exact decoding

Every datum of every covered format is exactly representable in both
`Float32` and `Float64`; `FloatBytes.decode` is a table lookup into a
precomputed atlas and is exact, never rounded (`Float64(x)` and
`Float32(x)` expose it). This is asserted
exhaustively by the test suite, not assumed.

```jldoctest formats
julia> Float64(FloatBytes.fromcode(Binary8p1uf, 0x01))   # 2^-127, exact
5.877471754111438e-39
```

## Format-level queries

Draft §4.14 queries, all folding to constants at compile time:
[`bitwidth`](@ref), `Base.precision`, [`signedness`](@ref),
[`domain`](@ref), [`exponent_bitwidth`](@ref),
[`trailing_significand_bits`](@ref), [`exponent_bias`](@ref); value
queries [`maxfinite`](@ref), [`minfinite`](@ref),
[`minpositive`](@ref), [`maxsubnormal`](@ref) (NaN when `P == 1` — no
subnormals exist), [`minnormal`](@ref).
