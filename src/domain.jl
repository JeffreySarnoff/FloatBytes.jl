# domain.jl — the Byte Domain (ImplementationPlan §5).
# FloatByte type, traits, code identity, format lattice, decode atlas.

abstract type Signedness end
struct SignedFormat   <: Signedness end
struct UnsignedFormat <: Signedness end
abstract type DomainKind end
struct FiniteDomain   <: DomainKind end
struct ExtendedDomain <: DomainKind end

"""
    FloatByte{K,P,S<:Signedness,D<:DomainKind} <: AbstractFloat

One-byte floating-point value per IEEE P3109 draft 1: bitwidth `K ∈ 3:8`,
precision `P`, signedness marker `S`, domain marker `D`. Users normally
use the 120 exported aliases (`Binary8p4se`, …) rather than the raw
parameters. Canonical-byte invariant: bits above `K` are zero for every
native value.
"""
primitive type FloatByte{K,P,S<:Signedness,D<:DomainKind} <: AbstractFloat 8 end

# --- traits (all fold at compile time; gate G9) ----------------------------

kbits(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = K
pbits(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = P
issigned(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = S === SignedFormat
isextended(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = D === ExtendedDomain

# CRITICAL: all shift/subtract arithmetic in Int, narrowed only at the end.
# `UInt8(1) << 8 == 0x00` would silently alias K = 8 codes (§1.2).
bias(::Type{T}) where {T<:FloatByte} =
    issigned(T) ? 1 << (kbits(T) - pbits(T) - 1) : 1 << (kbits(T) - pbits(T))
expbits(::Type{T}) where {T<:FloatByte} =
    issigned(T) ? kbits(T) - pbits(T) : kbits(T) - pbits(T) + 1
trailingbits(::Type{T}) where {T<:FloatByte} = pbits(T) - 1

# Special code points (§1.2, Annex B Table 3). Absent specials are `nothing`.
nan_code(::Type{T}) where {T<:FloatByte} =
    UInt8(issigned(T) ? 1 << (kbits(T) - 1) : (1 << kbits(T)) - 1)
pinf_code(::Type{T}) where {T<:FloatByte} =
    isextended(T) ? UInt8(issigned(T) ? (1 << (kbits(T) - 1)) - 1 : (1 << kbits(T)) - 2) : nothing
ninf_code(::Type{T}) where {T<:FloatByte} =
    (isextended(T) && issigned(T)) ? UInt8((1 << kbits(T)) - 1) : nothing
one_code(::Type{T}) where {T<:FloatByte} =
    UInt8(issigned(T) ? 1 << (kbits(T) - 2) : 1 << (kbits(T) - 1))
zero_code(::Type{T}) where {T<:FloatByte} = 0x00

# Maximum-finite code (magnitude part for signed).
function maxfinite_code(::Type{T}) where {T<:FloatByte}
    K = kbits(T)
    if issigned(T)
        return UInt8(isextended(T) ? (1 << (K - 1)) - 2 : (1 << (K - 1)) - 1)
    else
        return UInt8(isextended(T) ? (1 << K) - 3 : (1 << K) - 2)
    end
end
minpositive_code(::Type{T}) where {T<:FloatByte} = 0x01

"""
    validate_format(T)

Check `K ∈ 3:8` and the legal `P` range for the signedness (§1.1);
folds away entirely for constant formats.
"""
function validate_format(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D}
    ok = (K isa Int) && (P isa Int) && 3 <= K <= 8 &&
         (S === SignedFormat ? (1 <= P <= K - 1) : (1 <= P <= K))
    ok || throw(ArgumentError("invalid FloatByte parameters K=$K P=$P Σ=$S Δ=$D"))
    return nothing
end

# --- code identity (§5.3) --------------------------------------------------

"""
    fromcode(T, code::Integer) -> T

The only safe code-identity constructor. Range-checked; `reinterpret`
into a `FloatByte` is documented unsafe and confined to this file and the
dense execution lowering.
"""
function fromcode(::Type{T}, code::Integer) where {T<:FloatByte}
    validate_format(T)
    (0 <= code < (1 << kbits(T))) ||
        throw(ArgumentError("code $code out of range for $(T)"))
    return reinterpret(T, UInt8(code))
end

codepoint(x::FloatByte) = reinterpret(UInt8, x)
isvalidcode(::Type{T}, code::Integer) where {T<:FloatByte} =
    0 <= code < (1 << kbits(T))

# --- code classification (pure code-point predicates) ----------------------

const CODE_FINITE = 0x00
const CODE_ZERO   = 0x01   # finite as well; kept distinct for kernels
const CODE_NAN    = 0x02
const CODE_PINF   = 0x03
const CODE_NINF   = 0x04

@inline function classify_code(::Type{T}, c::UInt8) where {T<:FloatByte}
    c == 0x00 && return CODE_ZERO
    c == nan_code(T) && return CODE_NAN
    if isextended(T)
        c == pinf_code(T) && return CODE_PINF
        issigned(T) && c == ninf_code(T) && return CODE_NINF
    end
    return CODE_FINITE
end

@inline is_nan_code(::Type{T}, c::UInt8) where {T<:FloatByte} = c == nan_code(T)

# --- field decode (§1.3) ---------------------------------------------------

"""
    decode_fields(T, code) -> (S, E, neg)

Canonical `(significand, exponent, sign)` of a *finite* code:
value `= (-1)^neg * S * 2^(E + 1 - P)`, with `2^(P-1) ≤ S < 2^P` normal
and `0 ≤ S < 2^(P-1)` subnormal at `E = 1 - B`. Pure bit arithmetic; the
caller has already excluded special codes.
"""
@inline function decode_fields(::Type{T}, code::UInt8) where {T<:FloatByte}
    K = kbits(T); P = pbits(T); B = bias(T)
    c = Int(code)
    if issigned(T)
        neg = c >= (1 << (K - 1))
        mag = c & ((1 << (K - 1)) - 1)
    else
        neg = false
        mag = c
    end
    e = mag >> (P - 1)
    t = mag & ((1 << (P - 1)) - 1)
    if e == 0
        return (t, 1 - B, neg)                  # subnormal (S == 0 iff zero)
    else
        return (t + (1 << (P - 1)), e - B, neg) # normal
    end
end

"""
    unpack(T, code) -> Datum

Exact internal form of a finite nonzero code. Specials and zero are
excluded by the caller. Must not read the atlas (pure bit arithmetic).
"""
@inline function unpack(::Type{T}, code::UInt8) where {T<:FloatByte}
    S, E, neg = decode_fields(T, code)
    return Datum(UInt64(S), Int32(E + 1 - pbits(T)), neg)
end

# Maximum finite value as canonical (S, E); used by saturation.
function maxfinite_fields(::Type{T}) where {T<:FloatByte}
    S, E, _ = decode_fields(T, maxfinite_code(T))
    return (S, E)
end
emax(::Type{T}) where {T<:FloatByte} = maxfinite_fields(T)[2]

# --- format lattice (§1.5) -------------------------------------------------

# K ascending; Signed before Unsigned; P ascending; Finite before Extended.
# 4K - 2 formats per K; 120 total; 13 296 code points.
const FORMAT_ID_BASE = (0, 10, 24, 42, 64, 90)   # formats before each K ∈ 3:8

function format_id(K::Int, P::Int, signed::Bool, extended::Bool)
    base = FORMAT_ID_BASE[K - 2]
    within = signed ? (P - 1) * 2 + (extended ? 1 : 0) :
                      2 * (K - 1) + (P - 1) * 2 + (extended ? 1 : 0)
    return base + within + 1
end
format_id(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} =
    format_id(K, P, S === SignedFormat, D === ExtendedDomain)

"""
    format_params(id) -> (K, P, signed, extended)

Inverse of `format_id`.
"""
function format_params(id::Int)
    1 <= id <= 120 || throw(ArgumentError("format id $id out of 1:120"))
    K = 3
    while K < 8 && id > FORMAT_ID_BASE[K - 2] + 4 * K - 2
        K += 1
    end
    within = id - FORMAT_ID_BASE[K - 2] - 1
    signed = within < 2 * (K - 1)
    signed || (within -= 2 * (K - 1))
    P = (within >> 1) + 1
    extended = isodd(within)
    return (K, P, signed, extended)
end

function format_type(id::Int)
    K, P, s, e = format_params(id)
    return FloatByte{K, P, s ? SignedFormat : UnsignedFormat,
                     e ? ExtendedDomain : FiniteDomain}
end

# FORMAT_OFFSETS[id] = atlas offset in code points; last entry 13 296.
const FORMAT_OFFSETS = let
    offs = Memory{Int}(undef, 121)
    offs[1] = 0
    for id in 1:120
        K, _, _, _ = format_params(id)
        offs[id + 1] = offs[id] + (1 << K)
    end
    offs
end
const N_CODEPOINTS = FORMAT_OFFSETS[end]

# --- the decode atlas (§5.4) -----------------------------------------------

function _decode_f64(::Type{T}, code::UInt8) where {T<:FloatByte}
    cls = classify_code(T, code)
    cls == CODE_NAN  && return NaN
    cls == CODE_PINF && return Inf
    cls == CODE_NINF && return -Inf
    S, E, neg = decode_fields(T, code)
    v = S * exp2(E + 1 - pbits(T))              # exact: every datum is a
    return neg ? -v : v                         # Float64-representable dyadic
end

function _build_atlas()
    a = Memory{Float64}(undef, N_CODEPOINTS)
    for id in 1:120
        F = format_type(id)
        off = FORMAT_OFFSETS[id]
        for c in 0:(1 << kbits(F)) - 1
            a[off + c + 1] = _decode_f64(F, UInt8(c))
        end
    end
    return a
end

const DECODE64 = _build_atlas()
const DECODE32 = let
    a = Memory{Float32}(undef, N_CODEPOINTS)
    for i in 1:N_CODEPOINTS
        a[i] = Float32(DECODE64[i])             # exact: T1 asserts both atlases
    end
    a
end

"""
    decode(x::FloatByte) -> Float64

Exact decoded value (`NaN`, `±Inf` for specials). One const-Memory load
plus a bounds-elided indexed read; the offset folds per alias.
"""
@inline function decode(x::FloatByte{K,P,S,D}) where {K,P,S,D}
    T = FloatByte{K,P,S,D}
    @inbounds DECODE64[FORMAT_OFFSETS[format_id(T)] + Int(codepoint(x)) + 1]
end
@inline function decode32(x::FloatByte{K,P,S,D}) where {K,P,S,D}
    T = FloatByte{K,P,S,D}
    @inbounds DECODE32[FORMAT_OFFSETS[format_id(T)] + Int(codepoint(x)) + 1]
end
