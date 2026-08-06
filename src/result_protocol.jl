# result_protocol.jl — the tier-1 window substrate (ImplementationPlan §4).
# Loads first; imports nothing else from the package. Everything here is
# isbits and allocation-free. This is the *shipped* exact-evaluation path;
# the arbitrary-precision path lives in reference/ and is never reachable
# from src/.

const WINDOW_LSB = 40            # "G" in ComposedPlan §1.1

"""
    Datum

A finite nonspecial value `(-1)^neg * sig * 2^exp`, exactly.
Operand significands are < 2^8; products < 2^16; the four-factor dot
products of §13.1 stay < 2^32. `sig == 0` iff the value is zero.
`Datum` never holds a special value.
"""
struct Datum
    sig::UInt64
    exp::Int32
    neg::Bool
end

"""
    Aligned

A magnitude known to within one window unit:
`magnitude = (w + δ) * 2^exp` with `δ == 0` iff `tail == false`
and `δ ∈ (0, 1)` (strictly) iff `tail == true`.

The correctness contract with `project_code` (§7.1): when `tail` is set,
no rounding boundary of the target projection may lie strictly inside the
open interval `(w, w + 1) * 2^exp`. Every kernel in this file guarantees
that by carrying enough significant bits (≥ 33 fraction bits below the
leading bit — deterministic modes and stochastic `N ≤ 32`).
"""
struct Aligned
    w::UInt64
    exp::Int32
    neg::Bool
    tail::Bool
end

negate(d::Datum) = Datum(d.sig, d.exp, !d.neg)

@inline topexp(d::Datum) = Int(d.exp) + (63 - leading_zeros(d.sig))   # requires sig ≠ 0

# Order operands for the window anchor. Lexicographic (exp, sig) is exact on
# equal exponents (the wraparound-critical case); when exponents differ it can
# misorder close magnitudes, which the exact branch of `window_sum` tolerates
# by computing the symmetric difference.
@inline cmp_mag(x::Datum, y::Datum) = cmp((x.exp, x.sig), (y.exp, y.sig))

"""
    window_sum(x, y) -> Aligned

Two-term exact-or-tailed sum. Caller guarantees `x.sig ≠ 0 ≠ y.sig`
(zeros short-circuit in the kernel). See ImplementationPlan §4.2 for the
overflow and separation bounds:
  * `hi.sig ≤ 2^16 - 1` so `w ≤ 2^56` and `w + v < 2^57` — no overflow;
  * in the drop branch the separation exceeds 24 binades, so
    `v < w / 2^24` and `w - v - 1 > 0`.
"""
function window_sum(x::Datum, y::Datum)::Aligned
    hi, lo = cmp_mag(x, y) >= 0 ? (x, y) : (y, x)
    w  = hi.sig << WINDOW_LSB
    e0 = Int(hi.exp) - WINDOW_LSB
    s  = Int(lo.exp) - e0                      # lo's shift into the window; ≤ 40
    if s >= 0                                  # retained exactly
        v = lo.sig << s
        if hi.neg == lo.neg
            return Aligned(w + v, Int32(e0), hi.neg, false)
        elseif w >= v
            # equal-magnitude cancellation gives w == v: exact zero, sign +
            # is enforced at projection (w == 0 encodes code 0).
            return Aligned(w - v, Int32(e0), hi.neg, false)
        else
            # cmp_mag misordered close magnitudes across different exponents;
            # both terms are exact in this frame, so flip the difference.
            return Aligned(v - w, Int32(e0), !hi.neg, false)
        end
    else                                       # partial or total drop
        sh   = -s
        v    = sh >= 64 ? UInt64(0) : lo.sig >> sh
        tail = sh >= 64 ? lo.sig != 0 : (lo.sig & ((UInt64(1) << sh) - 1)) != 0
        if hi.neg == lo.neg
            return Aligned(w + v, Int32(e0), hi.neg, tail)
        elseif tail                            # w - v - ε ≡ (w - v - 1) + (1 - ε)
            @assert w > v + 1
            return Aligned(w - v - 1, Int32(e0), hi.neg, true)
        else
            return Aligned(w - v, Int32(e0), hi.neg, false)
        end
    end
end

@inline window_mul(x::Datum, y::Datum) =
    Datum(x.sig * y.sig, x.exp + y.exp, xor(x.neg, y.neg))   # ≤ 16 bits, exact

@inline window_fma(x::Datum, y::Datum, z::Datum) = window_sum(window_mul(x, y), z)

"""
    window_div(x, y) -> Aligned

`x / y` for finite nonzero datums. The 56-bit-scaled quotient carries at
least 49 significant bits; the residue `r / y.sig` is not dyadic in
general, but it lies strictly inside one quotient unit, so the sticky bit
is exact for every deterministic mode and stochastic `N ≤ 32`
(no rounding boundary can sit strictly inside an open unit interval —
boundaries are integer multiples of ≥ 2^16 units at this width).
"""
@inline function window_div(x::Datum, y::Datum)::Aligned
    q, r = divrem(x.sig << 56, y.sig)          # x.sig ≤ 255: no overflow
    Aligned(q, Int32(Int(x.exp) - 56 - Int(y.exp)), xor(x.neg, y.neg), r != 0)
end

# --- 128-bit integer square root -------------------------------------------

"""
    isqrt128(m::UInt128) -> UInt64

Floor square root, exact, for `m < 2^126`. Float64 seed, two Newton steps,
then bounded final corrections (each loop runs at most a couple of times).
"""
function isqrt128(m::UInt128)::UInt64
    m == 0 && return UInt64(0)
    r = UInt128(min(trunc(UInt64, sqrt(Float64(m))), typemax(UInt64) >> 1))
    r == 0 && (r = UInt128(1))
    r = (r + m ÷ r) >> 1
    r = (r + m ÷ r) >> 1
    while r * r > m
        r -= 1
    end
    while (r + 1) * (r + 1) <= m
        r += 1
    end
    return UInt64(r)
end

"""
    window_sqrt(x) -> Aligned

Square root of an exact datum; `x.neg == false`, `x.sig ≠ 0`.
Widened beyond ImplementationPlan §4.6 to a 2^84 scale so the result
carries ≥ 42 significant bits: the stochastic `N ≤ 32` predicates then
read exact top-`N+1` fraction bits directly from the window.
"""
function window_sqrt(x::Datum)::Aligned
    sh = iseven(Int(x.exp) - 84) ? 84 : 85     # make (exp - sh) even
    m  = UInt128(x.sig) << sh                  # ≤ 255·2^85 < 2^93
    s  = isqrt128(m)                           # ∈ [2^42, 2^47)
    Aligned(s, Int32((Int(x.exp) - sh) >> 1), false, widemul(s, s) != m)
end

"""
    window_rsqrt(x) -> Aligned

Direct reciprocal square root (never composed div∘sqrt: that would feed a
tailed divisor into a routine whose reasoning assumes exact operands).
Widened to a 2^94 target so the result carries ≥ 42 significant bits.
The at-most-one correction step in each direction is proved in
ImplementationPlan §4.7 (the proof accounts for the flooring `Q ÷ m`).
"""
function window_rsqrt(x::Datum)::Aligned
    m, e = iseven(Int(x.exp)) ? (x.sig, Int(x.exp)) : (x.sig << 1, Int(x.exp) - 1)
    Q  = UInt128(1) << 94
    w  = isqrt128(Q ÷ UInt128(m))              # candidate; w ∈ [2^42, 2^47]
    m128 = UInt128(m)
    if widemul(w + one(UInt64), w + one(UInt64)) * m128 <= Q
        w += one(UInt64)
    end
    if widemul(w, w) * m128 > Q
        w -= one(UInt64)
    end
    @assert widemul(w + one(UInt64), w + one(UInt64)) * m128 > Q
    @assert widemul(w, w) * m128 <= Q
    Aligned(w, Int32(-47 - (e >> 1)), false, widemul(w, w) * m128 != Q)
end

"""
    window_hypot(x, y) -> Aligned

`sqrt(x^2 + y^2)` for finite nonzero datums. The squares are exact 16-bit
products, both non-negative — `window_sum` applies with no cancellation.
The root is then taken at ~2^116 scale (≥ 49 significant bits); the true
root lies strictly inside one unit of the result frame unless it is exact,
so the sticky representative decides every deterministic mode and every
stochastic mode with `N ≤ 32` (ImplementationPlan §4.8, simplified: at
this width no rounding boundary can be interior to the open unit
interval, so the step-5 boundary refinement is never needed).
"""
function window_hypot(x::Datum, y::Datum)::Aligned
    a = window_sum(window_mul(x, x), window_mul(y, y))
    a.w == 0 && return Aligned(UInt64(0), Int32(0), false, false)
    sh = iseven(Int(a.exp) - 58) ? 58 : 59
    W  = UInt128(a.w) << sh                    # ≤ 2^57·2^59 = 2^116
    s  = isqrt128(W)                           # ≥ 2^49
    e2 = Int32((Int(a.exp) - sh) >> 1)
    exact = widemul(s, s) == W && !a.tail
    Aligned(s, e2, false, !exact)
end

# --- fixed-limb signed superaccumulator (ImplementationPlan §13.1) ---------
#
# Shipped here rather than in blocks.jl because the FAA kernel (§4.4) is the
# B = 3 specialization of this machinery and result_protocol must stay
# self-contained. Exact by construction: integer addition into an L-limb
# fixed-point frame captures every bit.

struct LimbAcc{L}
    pos::NTuple{L, UInt64}   # non-negative running sum, fixed-point at 2^emin
    neg::NTuple{L, UInt64}   # non-negative subtractor; value = pos - neg
end

@inline ntuple_zero(::Val{L}) where {L} = ntuple(_ -> UInt64(0), Val(L))

LimbAcc{L}() where {L} = LimbAcc{L}(ntuple_zero(Val(L)), ntuple_zero(Val(L)))

# Carry-propagating limb addition — recursive tail-call form; unrolled at
# compile time because L is a type parameter. The naive closure-over-carry
# form boxes the carry (captured-variable rules) and is forbidden.
@inline add_limbs(a::NTuple{L, UInt64}, b::NTuple{L, UInt64}) where {L} =
    _add_limbs(a, b, UInt64(0), Val(1), Val(L))

@inline function _add_limbs(a, b, carry::UInt64, ::Val{I}, ::Val{L}) where {I, L}
    s,  c1 = Base.add_with_overflow(a[I], b[I])
    s2, c2 = Base.add_with_overflow(s, carry)
    ncarry = UInt64(c1) + UInt64(c2)
    if I == L
        return (s2,)
    else
        return (s2, _add_limbs(a, b, ncarry, Val(I + 1), Val(L))...)
    end
end

@inline sub_limbs(a::NTuple{L, UInt64}, b::NTuple{L, UInt64}) where {L} =
    _sub_limbs(a, b, UInt64(0), Val(1), Val(L))

@inline function _sub_limbs(a, b, borrow::UInt64, ::Val{I}, ::Val{L}) where {I, L}
    d,  b1 = Base.sub_with_overflow(a[I], b[I])
    d2, b2 = Base.sub_with_overflow(d, borrow)
    nborrow = UInt64(b1) + UInt64(b2)
    if I == L
        return (d2,)
    else
        return (d2, _sub_limbs(a, b, nborrow, Val(I + 1), Val(L))...)
    end
end

# Limb comparison, most-significant limb first.
@inline function cmp_limbs(a::NTuple{L, UInt64}, b::NTuple{L, UInt64}) where {L}
    i = L
    while i >= 1
        a[i] != b[i] && return a[i] < b[i] ? -1 : 1
        i -= 1
    end
    return 0
end

# Sparse addend: sig placed at bit offset `off` in the L-limb frame.
@inline function sparse_addend(::Val{L}, off::Int, sig::UInt64) where {L}
    limb  = (off >> 6) + 1                     # 1-based limb of the LSB
    intra = off & 0x3f
    lo_word = sig << intra
    hi_word = intra == 0 ? UInt64(0) : sig >> (64 - intra)
    ntuple(i -> i == limb ? lo_word : (i == limb + 1 ? hi_word : UInt64(0)), Val(L))
end

"""
    update(acc, d, emin) -> LimbAcc

Add exact `Datum` `d` into the accumulator whose fixed-point frame has its
LSB at `2^emin`. `d.exp - emin ≥ 0` by construction of the frame.
"""
@inline function update(acc::LimbAcc{L}, d::Datum, emin::Int) where {L}
    d.sig == 0 && return acc
    off = Int(d.exp) - emin
    @assert off >= 0
    addend = sparse_addend(Val(L), off, d.sig)
    if d.neg
        return LimbAcc{L}(acc.pos, add_limbs(acc.neg, addend))
    else
        return LimbAcc{L}(add_limbs(acc.pos, addend), acc.neg)
    end
end

# Canonicalizing merge: exactly one of {pos, neg} nonzero afterwards.
@inline function merge_acc(a::LimbAcc{L}, b::LimbAcc{L}) where {L}
    pos = add_limbs(a.pos, b.pos)
    neg = add_limbs(a.neg, b.neg)
    return cmp_limbs(pos, neg) >= 0 ?
        LimbAcc{L}(sub_limbs(pos, neg), ntuple_zero(Val(L))) :
        LimbAcc{L}(ntuple_zero(Val(L)), sub_limbs(neg, pos))
end

"""
    close_to_aligned(acc, emin) -> Aligned

Resolve `|pos - neg|` and collapse to the top-64-bit window with a sticky
OR of every dropped lower bit. The window then carries `64 - P ≥ 56`
exact fraction bits below its leading bit — comfortably above the 33-bit
requirement of every deterministic mode and stochastic `N ≤ 32`.
"""
function close_to_aligned(acc::LimbAcc{L}, emin::Int)::Aligned where {L}
    c = cmp_limbs(acc.pos, acc.neg)
    c == 0 && return Aligned(UInt64(0), Int32(0), false, false)
    neg = c < 0
    mag = neg ? sub_limbs(acc.neg, acc.pos) : sub_limbs(acc.pos, acc.neg)
    j = L
    while mag[j] == 0
        j -= 1
    end
    h = 63 - leading_zeros(mag[j])
    w = mag[j] << (63 - h)
    tail = false
    if j >= 2
        w |= mag[j - 1] >> (h + 1)                    # h + 1 ∈ 1:64; >> 64 == 0
        tail = (mag[j - 1] << (63 - h)) != 0          # the h + 1 dropped low bits
        i = j - 2
        while i >= 1 && !tail
            tail = mag[i] != 0
            i -= 1
        end
    end
    exp = emin + (j - 1) * 64 + h - 63
    return Aligned(w, Int32(exp), neg, tail)
end
