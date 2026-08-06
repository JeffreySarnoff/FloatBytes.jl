# external.jl — Float64/Float32 bridges (§14.3).
# Bridge semantics (draft §4.8): external decode collapses every host NaN
# to the single internal NaN and both host zero signs to the single zero;
# external encode emits positive zero and a quiet NaN with zero payload.
# The bit-level path below is exact for every finite host value: single
# projection, no double rounding.

# Exact (mantissa, exponent, sign) of a finite nonzero host float.
@inline function _float_fields(x::Float64)
    u = reinterpret(UInt64, x)
    neg = (u >> 63) != 0
    ebits = Int((u >> 52) & 0x7ff)
    frac = u & 0x000fffffffffffff
    if ebits == 0
        return (frac, -1074, neg)                       # subnormal
    else
        return (frac | (UInt64(1) << 52), ebits - 1075, neg)
    end
end
@inline function _float_fields(x::Float32)
    u = reinterpret(UInt32, x)
    neg = (u >> 31) != 0
    ebits = Int((u >> 23) & 0xff)
    frac = UInt64(u & 0x007fffff)
    if ebits == 0
        return (frac, -149, neg)
    else
        return (frac | (UInt64(1) << 23), ebits - 150, neg)
    end
end

"""
    quantize(F, x, ρ = DEFAULT_PROJECTION; entropy = nothing, index = 0) -> F

Mathematical construction: project real `x` into format `F` under ρ.
This is the visible projecting constructor; `fromcode` is the raw one.
"""
function quantize(::Type{F}, x::Union{Float64, Float32},
                  ρ::Projection = DEFAULT_PROJECTION;
                  entropy = nothing, index::Integer = 0) where {F<:FloatByte}
    validate_format(F)
    e = _entropy_word_for(ρ, entropy, index)
    isnan(x) && return fromcode(F, nan_code(F))
    if isinf(x)
        return fromcode(F, project_special(F, ρ, x > 0 ? CODE_PINF : CODE_NINF))
    end
    iszero(x) && return fromcode(F, 0x00)               # −0.0 collapses
    m, ex, neg = _float_fields(x)
    return fromcode(F, project_code(F, ρ, Aligned(m, Int32(ex), neg, false), e))
end

function quantize(::Type{F}, x::Real, ρ::Projection = DEFAULT_PROJECTION;
                  entropy = nothing, index::Integer = 0) where {F<:FloatByte}
    validate_format(F)
    e = _entropy_word_for(ρ, entropy, index)
    if x isa Integer
        iszero(x) && return fromcode(F, 0x00)
        neg = x < 0
        b = BigInt(abs(x))
        return fromcode(F, _project_bigint_scaled(F, ρ, b, 0, neg, e))
    end
    # Exact slow path via Rational{BigInt}; correctness over speed for
    # uncommon input types (Rational, BigFloat, Float16, ...).
    if x isa BigFloat
        isnan(x) && return fromcode(F, nan_code(F))
        isinf(x) && return fromcode(F, project_special(F, ρ, x > 0 ? CODE_PINF : CODE_NINF))
        iszero(x) && return fromcode(F, 0x00)
    end
    r = rationalize_exact(x)
    iszero(r) && return fromcode(F, 0x00)
    neg = r < 0
    r = abs(r)
    # Scale numerator to a 63..65-bit quotient, take the 64-bit window and
    # a sticky for the remainder: single projection.
    n = numerator(r); d = denominator(r)
    sh = 64 + ndigits(d, base = 2) - ndigits(n, base = 2)
    _scaled_divrem(n, d, sh) = sh >= 0 ? divrem(n << sh, d) : divrem(n, d << (-sh))
    q, rem = _scaled_divrem(n, d, sh)
    while q >= (BigInt(1) << 64)
        sh -= 1
        q, rem = _scaled_divrem(n, d, sh)
    end
    return fromcode(F, project_code(F, ρ,
        Aligned(UInt64(q), Int32(-sh), neg, rem != 0), e))
end

rationalize_exact(x::AbstractFloat) = Rational{BigInt}(x)
rationalize_exact(x::Rational) = Rational{BigInt}(x)
rationalize_exact(x::Integer) = Rational{BigInt}(x)

function _project_bigint_scaled(::Type{F}, ρ::Projection, b::BigInt,
                                ex::Int, neg::Bool, e::UInt64) where {F<:FloatByte}
    nb = ndigits(b, base = 2)
    if nb <= 64
        return project_code(F, ρ, Aligned(UInt64(b), Int32(ex), neg, false), e)
    end
    sh = nb - 64
    w = UInt64(b >> sh)
    tail = b != (BigInt(w) << sh)
    return project_code(F, ρ, Aligned(w, Int32(ex + sh), neg, tail), e)
end

# --- external encode (exact by the atlas exactness theorem) ----------------

Base.Float64(x::FloatByte) = decode(x)
Base.Float32(x::FloatByte) = decode32(x)
Base.Float16(x::FloatByte) = Float16(decode(x))         # may round: Float16 is narrower
Base.BigFloat(x::FloatByte) = BigFloat(decode(x))       # exact

(::Type{F})(x::Real) where {F<:FloatByte} = quantize(F, x)
