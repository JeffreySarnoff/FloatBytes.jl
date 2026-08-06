```@meta
CurrentModule = FloatBytes
DocTestSetup = :(using FloatBytes)
```

# Julia interoperability

`FloatByte <: AbstractFloat`, with a deliberately curated surface: what
is defined is tested against the draft; what is not defined errors
loudly instead of silently borrowing host-float semantics.

## Arithmetic and promotion

Base operators (`+ - * / sqrt inv fma hypot abs min max clamp copysign`
and the transcendental names) work on **same-format** operands under the
scoped default projection. Cross-format arithmetic through Base
operators is an error pointing at the explicit seam; mixing with host
floats promotes to `Float64`/`Float32` (exact for every datum):

```jldoctest interop
julia> using FloatBytes

julia> Binary8p4se(1.5) + 1.0          # promotes to Float64, exactly
2.5

julia> Binary8p4se(1.5) + Binary8p3se(1.0)
ERROR: ArgumentError: cross-format + between Binary8p4se and Binary8p3se is not defined; use apply(Add, Fr, ρ, x, y) or Convert
[...]
```

## Host-float bridges

Decode is exact in both directions of the bridge; every host NaN
collapses to the single NaN, both host zeros to the single zero, and
construction from a host float is a single projection (no double
rounding):

```jldoctest interop
julia> Binary8p4se(-0.0) === zero(Binary8p4se)
true

julia> Float64(maxfinite(Binary8p4se))
224.0

julia> quantize(Binary8p4se, big"1.00000000000000000001")
Binary8p4se(0x40 ↦ 1.0)
```

## Ordering, sorting, hashing

`isless` follows the P3109 total order: **NaN sorts first** — one order
for predicates, sorting, and selection:

```jldoctest interop
julia> v = [Binary8p4se(2.0), Binary8p4se(NaN), Binary8p4se(-1.0)];

julia> sort(v)
3-element Vector{Binary8p4se}:
 Binary8p4se(0x80 ↦ NaN)
 Binary8p4se(0xc0 ↦ -1.0)
 Binary8p4se(0x48 ↦ 2.0)
```

`==`, `<`, `<=` follow the draft comparison predicates (NaN compares
false with everything, including itself); `hash` is on
`(format, code point)` so `Dict` and `Set` behave.

## Random values

`rand(rng, F)` draws uniformly over the **finite datums** of the format
(a documented choice — never NaN or `±Inf`).

## Zero-copy bytes

Arrays are plain byte buffers: `reinterpret(UInt8, ::Vector{F})` is
zero-copy, `sizeof(F) == 1`, `isbitstype(F)`.
