```@meta
CurrentModule = FloatBytes
```

# FloatBytes.jl

Byte-sized floating point for Julia 1.12+, implementing **IEEE P3109
working draft D1**: all 120 `Binary{K}p{P}{s|u}{e|f}` formats for
bitwidths `K ∈ 3:8`, with exact integer-window arithmetic, draft-literal
special-value semantics, every draft rounding mode (including three
stochastic variants), and saturation as a first-class component of the
projection specification.

!!! warning "Draft standard"
    P3109 is a working **draft**. This package tracks the repository's
    pinned draft revision (`FloatBytes.DRAFT_REVISION`) and makes **no
    conformance claim**. Semantics can change as the draft evolves.

## Installation

```julia
pkg> add https://github.com/JeffreySarnoff/FloatBytes.jl
```

Requires Julia 1.12 or later.

## Quick start

```jldoctest quick
julia> using FloatBytes

julia> x = Binary8p4se(1.5)              # projecting constructor (default ρ)
Binary8p4se(0x44 ↦ 1.5)

julia> x + Binary8p4se(2.25)             # Base operators, same format only
Binary8p4se(0x4f ↦ 3.75)

julia> apply(Divide, Binary8p4se, Projection(TowardZero(), SatFinite()), x, x + x)
Binary8p4se(0x38 ↦ 0.5)

julia> one(Binary8p4se) / zero(Binary8p4se)   # draft semantics: x/0 = NaN
Binary8p4se(0x80 ↦ NaN)

julia> Float64(x)                        # every datum is exact in Float64
1.5
```

Every value is one byte, `isbits`, and zero-copy reinterpretable to
`UInt8`:

```jldoctest quick
julia> a = Binary8p4se[1.0, 2.0, 3.0];

julia> reinterpret(UInt8, a)
3-element reinterpret(UInt8, ::Vector{Binary8p4se}):
 0x40
 0x48
 0x4c
```

## Where things differ from IEEE-754

All of these are normative in the draft, not implementation choices:

- **One NaN and one zero per format.** No NaN payloads, no `−0`. Host
  `−0.0` and every host NaN collapse at the bridge.
- **`x / 0 = NaN`** (not `±Inf`), and `Recip(±Inf) = 0`.
- **`RSqrt(0) = NaN`** (the draft defines `X ≤ 0 → NaN`).
- **NaN sorts first** in the total order — the opposite end from Float64
  sorting.
- **`NextGreaterThan(+Inf) = NaN`**, unlike IEEE-754 `nextUp`.
- **Finite-domain formats** and unsigned-negative results saturate per
  the projection's saturation mode (`SatNone` gives NaN where the draft
  says NaN).
- **Rounding precedes saturation**: values slightly above the maximum
  finite value round back into range before saturation classifies.

## Manual

```@contents
Pages = ["formats.md", "projection.md", "operations.md",
         "stochastic.md", "interop.md", "reference.md"]
Depth = 1
```

## Design documents

The implementation follows `docs/design/ImplementationPlan.md` in the
repository; deviations and stage coverage are recorded in
`docs/design/ImplementationStatus.md`. The independent verification
oracle lives in `reference/` and never shares rounding code with the
production path.
