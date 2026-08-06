# reference/draft_model.jl — independent draft transcription (§15).
# Written directly from IEEE_D1.md clauses; Rational{BigInt} values plus
# symbolic specials, plus an exact square-root value kind for the
# sqrt-family operations (⌊√r · 2^k⌋ = isqrt(⌊r · 4^k⌋), exact).
#
# ZERO imports from src/ (grep-enforced). Shares nothing with the
# production rounding path: production rounds through a 64-bit window with
# guard/round/sticky bits; this file decides the same questions from exact
# rational quotient/remainder arithmetic (§7.2 oracle discipline).

module DraftReference

export RFormat, all_formats, rcodes, decode_ref, project_ref, encode_ref,
       omega_add, omega_subtract, omega_multiply, omega_divide,
       omega_fma, omega_faa, omega_sqrt, omega_recip, omega_rsqrt,
       omega_hypot, RVal, SQ

struct RFormat
    K::Int
    P::Int
    signed::Bool
    extended::Bool
end

function all_formats()
    fs = RFormat[]
    for K in 3:8, signed in (true, false), P in 1:(signed ? K - 1 : K),
        extended in (false, true)
        push!(fs, RFormat(K, P, signed, extended))
    end
    return fs
end

rcodes(f::RFormat) = 0:(1 << f.K) - 1

rbias(f::RFormat) = f.signed ? 1 << (f.K - f.P - 1) : 1 << (f.K - f.P)
rnan(f::RFormat) = f.signed ? 1 << (f.K - 1) : (1 << f.K) - 1
rpinf(f::RFormat) = f.extended ? (f.signed ? (1 << (f.K - 1)) - 1 : (1 << f.K) - 2) : nothing
rninf(f::RFormat) = (f.extended && f.signed) ? (1 << f.K) - 1 : nothing
function rmaxfinite_code(f::RFormat)
    f.signed ? (f.extended ? (1 << (f.K - 1)) - 2 : (1 << (f.K - 1)) - 1) :
               (f.extended ? (1 << f.K) - 3 : (1 << f.K) - 2)
end

const Q = Rational{BigInt}

"""ωDecode per draft §4.7.2 → :nan | :pinf | :ninf | Rational{BigInt}"""
function decode_ref(f::RFormat, code::Integer)
    x = Int(code)
    x == rnan(f) && return :nan
    if f.extended
        x == rpinf(f) && return :pinf
        f.signed && x == rninf(f) && return :ninf
    end
    if f.signed && x >= (1 << (f.K - 1))
        v = decode_ref(f, x - (1 << (f.K - 1)))
        return -(v::Q)
    end
    T = BigInt(x) % (1 << (f.P - 1))
    Eb = BigInt(x) ÷ (1 << (f.P - 1))
    B = rbias(f)
    if Eb == 0
        return (0 + T // (BigInt(1) << (f.P - 1))) * Q(2)^(1 - B)
    else
        return (1 + T // (BigInt(1) << (f.P - 1))) * Q(2)^(Int(Eb) - B)
    end
end

# --- exact value kinds ------------------------------------------------------

"""SQ(r): the non-negative value √r for exact rational r ≥ 0."""
struct SQ
    r::Q
end

"""RVal: signed exact value — (neg, magnitude) where magnitude is Q or SQ."""
struct RVal
    neg::Bool
    mag::Union{Q, SQ}
end
RVal(v::Q) = RVal(v < 0, abs(v))

_is_zero(m::Q) = iszero(m)
_is_zero(m::SQ) = iszero(m.r)

# ⌊m · 2^k⌋ for m ≥ 0, any integer k. Exact.
function floor_scaled(m::Q, k::Int)::BigInt
    n, d = numerator(m), denominator(m)
    k >= 0 ? fld(n << k, d) : fld(n, d << (-k))
end
function floor_scaled(m::SQ, k::Int)::BigInt
    # ⌊√r · 2^k⌋ = ⌊√(r · 4^k)⌋ = isqrt(⌊r · 4^k⌋)  (proof: monotone squares)
    return isqrt(floor_scaled(m.r, 2k))
end

# Is m · 2^k an exact integer?
function is_int_scaled(m::Q, k::Int)::Bool
    n, d = numerator(m), denominator(m)
    k >= 0 ? (n << k) % d == 0 : n % (d << (-k)) == 0
end
function is_int_scaled(m::SQ, k::Int)::Bool
    r = m.r * Q(BigInt(1) << (2 * max(k, 0)), BigInt(1) << (2 * max(-k, 0)))
    denominator(r) == 1 || return false
    s = isqrt(numerator(r))
    return s * s == numerator(r)
end

# floor(log2 m) for m > 0: integer estimate, then exact fix-up so that
# 2^e ≤ m < 2^(e+1) (compares are exact via floor_scaled).
_est_log2(m::Q) = ndigits(numerator(m), base = 2) - ndigits(denominator(m), base = 2)
_est_log2(m::SQ) = fld(_est_log2(m.r), 2)
function floor_log2(m)::Int
    e = _est_log2(m)
    while floor_scaled(m, -e) < 1
        e -= 1
    end
    while floor_scaled(m, -(e + 1)) >= 1
        e += 1
    end
    return e
end

include("rational_projection.jl")

end # module
