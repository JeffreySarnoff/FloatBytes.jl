# projection.jl — the only code-producing implementation (§7).
# Exactly round_to_precision, then saturate, then encode — never merged:
# rounding precedes saturation, the rounding mode stays visible to
# saturation, and saturation never rounds (draft §4.7.3 NOTE 1, §4.7.5).

# --- rounding predicates (draft §4.7.4; ImplementationPlan §7.1) -----------
#
# Uniform inputs: `keep` (truncated significand ⌊S̃⌋), `rbit` (fraction bit
# ν ≥ 1/2), `sticky` (any fraction below rbit, including the tail), `neg`,
# and for stochastic modes the exact top-M fraction bits extracted by the
# caller. All pure integer arithmetic; ν is never converted to a float.

# CodeIsEven (draft): P > 1 → IsEven(⌊S̃⌋); P == 1 → ⌊S̃⌋ == 0 or IsEven(Q+B).
@inline _code_is_even(keep::UInt64, ulpexp::Int, P::Int, B::Int) =
    P > 1 ? iseven(keep) : (keep == 0 || iseven(ulpexp + B))

@inline round_away(::TowardZero, keep, rbit, sticky, neg, ulpexp, P, B) = false
@inline round_away(::TowardPositive, keep, rbit, sticky, neg, ulpexp, P, B) =
    !neg && (rbit | sticky)
@inline round_away(::TowardNegative, keep, rbit, sticky, neg, ulpexp, P, B) =
    neg && (rbit | sticky)
@inline round_away(::NearestTiesToAway, keep, rbit, sticky, neg, ulpexp, P, B) = rbit
@inline round_away(::NearestTiesToEven, keep, rbit, sticky, neg, ulpexp, P, B) =
    rbit && (sticky || !_code_is_even(keep, ulpexp, P, B))
@inline round_away(::ToOdd, keep, rbit, sticky, neg, ulpexp, P, B) =
    (rbit | sticky) && _code_is_even(keep, ulpexp, P, B)

# Exact floor(ν · 2^M) from the window: the k fraction bits of w below the
# ulp are exact, and the tail contributes strictly less than one unit at
# bit position k, so the floor is exact whenever M ≤ k (asserted at the
# call site). When the carrier is exact (no tail) the bits beyond k are
# genuine zeros and any M is fine.
@inline function _nu_scaled(w::UInt64, k::Int, M::Int)
    k <= 0 && return UInt64(0)
    low = k >= 64 ? w : w & ((UInt64(1) << k) - 1)
    M >= k && return low << (M - k)   # exact only for exact carriers (asserted)
    sh = k - M
    return sh >= 64 ? UInt64(0) : low >> sh
end

# entropy: raw 64-bit word; R = top N bits (0 ≤ R < 2^N).
@inline function round_away_stochastic(::StochasticA{N}, w, k, tail, entropy) where {N}
    t = _nu_scaled(w, k, N)
    R = entropy >> (64 - N)
    return t + R >= (UInt64(1) << N)
end
@inline function round_away_stochastic(::StochasticB{N}, w, k, tail, entropy) where {N}
    t = _nu_scaled(w, k, N + 1)
    R = entropy >> (64 - N)
    return t + (2 * R + 1) >= (UInt64(1) << (N + 1))
end
@inline function round_away_stochastic(::StochasticC{N}, w, k, tail, entropy) where {N}
    t = _nu_scaled(w, k, N)                    # integer part of ν·2^N
    # RNITE tie information: next bit and everything below it (plus tail).
    b = _nu_scaled(w, k, N + 1) & 0x01
    rest = tail
    if !rest && k > N + 1
        kk = k - N - 1
        mask = kk >= 64 ? typemax(UInt64) : (UInt64(1) << kk) - 1
        rest = (w & mask) != 0
    end
    r = t + ((b == 1 && (rest || isodd(t))) ? 1 : 0)
    R = entropy >> (64 - N)
    return r + R >= (UInt64(1) << N)
end

# --- round to precision (§7.1) ---------------------------------------------

"""
    round_to_precision(a::Aligned, P, B, mode, entropy)
        -> (keep::UInt64, E::Int, neg::Bool)   [keep == 0 encodes zero]

Round to `P` significant bits with unbounded exponent (subnormal quantum
clamp applied; overflow left in place for saturation to classify).
"""
@inline function round_to_precision(a::Aligned, P::Int, B::Int,
                                    mode::RoundingMode, entropy::UInt64)
    w = a.w
    h = 63 - leading_zeros(w)                  # w == 0 handled before this
    msbexp = Int(a.exp) + h
    ulpexp = max(msbexp, 1 - B) + 1 - P        # subnormal quantum clamp
    k = ulpexp - Int(a.exp)                    # ulp bit position within w
    local keep::UInt64, rbit::Bool, sticky::Bool
    if k <= 0
        keep = w << (-k)                       # exact multiple of the ulp
        rbit = false
        sticky = a.tail
    elseif k > 64
        keep = UInt64(0)
        rbit = false
        sticky = w != 0 || a.tail
    elseif k == 64
        keep = UInt64(0)
        rbit = (w >> 63) != 0
        sticky = (w & (typemax(UInt64) >> 1)) != 0 || a.tail
    else
        keep = w >> k
        rbit = ((w >> (k - 1)) & 0x01) != 0
        sticky = (w & ((UInt64(1) << (k - 1)) - 1)) != 0 || a.tail
    end
    away = if mode isa StochasticA || mode isa StochasticB || mode isa StochasticC
        N = _stochastic_bits(typeof(mode))
        # The window guarantee: tailed carriers arrive with ≥ 33 exact
        # fraction bits, so the top-(N+1) extraction is exact for N ≤ 32.
        @assert !a.tail || k >= N + 1 "stochastic N=$N exceeds carrier width k=$k"
        round_away_stochastic(mode, w, k, a.tail, entropy)
    else
        round_away(mode, keep, rbit, sticky, a.neg, ulpexp, P, B)
    end
    if away
        keep += 1
        if keep == UInt64(1) << P              # carry out of the precision
            keep >>= 1
            ulpexp += 1
        end
    end
    return (keep, ulpexp + P - 1, a.neg)       # canonical (S, E)
end

# --- saturate ∘ encode (§7.3, §7.4, §1.4) ----------------------------------

# Encode a canonical in-range (S, E, neg); §1.4. Strict precondition: the
# datum is already in the target's datum set — encode never rounds and
# never saturates.
@inline function encode_fields(::Type{T}, S::Int, E::Int, neg::Bool) where {T<:FloatByte}
    P = pbits(T); B = bias(T)
    local mag::Int
    if S == 0
        return 0x00
    elseif S < (1 << (P - 1))
        @assert E == 1 - B
        mag = S
    else
        mag = ((E + B) << (P - 1)) | (S - (1 << (P - 1)))
    end
    code = neg ? mag | (1 << (kbits(T) - 1)) : mag
    return UInt8(code)
end

@inline maxfinite_code_signed(::Type{T}, neg::Bool) where {T<:FloatByte} =
    neg ? maxfinite_code(T) | UInt8(1 << (kbits(T) - 1)) : maxfinite_code(T)

# Saturation decision table (draft §4.7.5, implemented literally, ordered).
# Input: rounded finite (keep, E, neg) with unbounded E. Infinite and NaN
# auxiliary results never reach here — see `project_special`.
@inline function saturate_encode(::Type{T}, ρ::Projection{R,SM},
                                 keep::UInt64, E::Int, neg::Bool) where {T<:FloatByte, R, SM}
    keep == 0 && return 0x00
    S = Int(keep)
    Smax, Emax = maxfinite_fields(T)
    over = E > Emax || (E == Emax && S > Smax)
    neg_unrepresentable = neg && !issigned(T)          # X < Mlo = 0 in unsigned
    if !over && !neg_unrepresentable
        return encode_fields(T, S, E, neg)             # Mlo ≤ X ≤ Mhi
    end
    # Below, X > Mhi (over && !neg), or X < Mlo (signed: over && neg;
    # unsigned: any negative).
    above = !neg
    if SM === SatFinite
        above && return maxfinite_code(T)
        return issigned(T) ? maxfinite_code_signed(T, true) : 0x00   # Mlo
    elseif SM === SatPropagate
        above && return maxfinite_code(T)
        return issigned(T) ? maxfinite_code_signed(T, true) : 0x00
    else # SatNone
        # Directed rows first: rounding direction requires a finite result.
        if above && (R === TowardZero || R === TowardNegative)
            return maxfinite_code(T)
        end
        if !above && (R === TowardZero || R === TowardPositive)
            return issigned(T) ? maxfinite_code_signed(T, true) : 0x00
        end
        if isextended(T)
            above && return pinf_code(T)::UInt8
            issigned(T) && return ninf_code(T)::UInt8
            return nan_code(T)                          # negative, unsigned extended
        end
        return nan_code(T)                              # finite domain
    end
end

"""
    project_special(T, ρ, kind) -> UInt8

Projection of a NaN or infinite auxiliary result (draft §4.7.3 + §4.7.5
rows on `±∞`/`NaN`).
"""
function project_special(::Type{T}, ρ::Projection{R,SM}, kind::UInt8) where {T<:FloatByte, R, SM}
    kind == CODE_NAN && return nan_code(T)
    if kind == CODE_PINF
        SM === SatFinite && return maxfinite_code(T)
        if SM === SatPropagate
            return isextended(T) ? pinf_code(T)::UInt8 : maxfinite_code(T)
        end
        return isextended(T) ? pinf_code(T)::UInt8 : nan_code(T)   # SatNone
    end
    kind == CODE_NINF || throw(ArgumentError("not a special kind: $kind"))
    if SM === SatFinite
        return issigned(T) ? maxfinite_code_signed(T, true) : 0x00
    elseif SM === SatPropagate
        (issigned(T) && isextended(T)) && return ninf_code(T)::UInt8
        return issigned(T) ? maxfinite_code_signed(T, true) : 0x00
    else # SatNone
        (issigned(T) && isextended(T)) && return ninf_code(T)::UInt8
        return nan_code(T)
    end
end

"""
    project_code(T, ρ, a::Aligned, entropy) -> UInt8

The hot pointwise path: round to precision, saturate, encode.
"""
@inline function project_code(::Type{T}, ρ::Projection{R,SM},
                              a::Aligned, entropy::UInt64) where {T<:FloatByte, R, SM}
    if a.w == 0
        # `w == 0 && tail` cannot reach projection from any tier-1 kernel.
        @assert !a.tail
        return 0x00
    end
    keep, E, neg = round_to_precision(a, pbits(T), bias(T), R(), entropy)
    return saturate_encode(T, ρ, keep, E, neg)
end

# Projection of an exact datum (Convert, selection results, bridges).
@inline project_datum(::Type{T}, ρ::Projection, d::Datum, entropy::UInt64) where {T<:FloatByte} =
    d.sig == 0 ? 0x00 :
    project_code(T, ρ, Aligned(d.sig, d.exp, d.neg, false), entropy)
