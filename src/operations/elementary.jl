# operations/elementary.jl — the 25 unary transcendentals plus ArcTan2 /
# ArcTan2Pi (§9). Draft special rows first (symbolic, code-point tests);
# finite arguments go through a Ziv-style rigorous enclosure loop with
# certified directed MPFR endpoints. MPFR appears only in this file and in
# reference/ (grep-enforced, §17.6). Hypot is algebraic and NOT here.
#
# Enclosure decision: with `lo ≤ f(x) ≤ hi` (strict on both sides unless
# exact), scale so the window has ~61 significant bits and take integer
# floors L, H. If L == H, the true value lies strictly inside one open
# unit interval of that frame; every rounding boundary (datum, midpoint,
# stochastic granule multiple) is an integer multiple of ≥ 2^20 frame
# units, so no boundary is interior and the sticky representative decides
# every mode (deterministic and stochastic N ≤ 32). Otherwise double the
# precision; at the ceiling throw InconclusiveProjectionError — never a
# midpoint, never a guess.

const ELEM_P_START   = 96
const ELEM_P_CEILING = 4096

function _project_exact_big(::Type{Fr}, ρ::Projection, e::UInt64,
                            v::BigFloat) where {Fr<:FloatByte}
    iszero(v) && return 0x00
    neg = v < 0
    av = abs(v)
    ex = exponent(av)
    s = 63 - ex
    t = ldexp(av, s)
    w = trunc(BigInt, t)
    @assert t == w "exact transcendental value with oversized mantissa"
    return project_code(Fr, ρ, Aligned(UInt64(w), Int32(-s), neg, false), e)
end

function _try_window(::Type{Fr}, ρ::Projection, e::UInt64,
                     lo::BigFloat, hi::BigFloat) where {Fr<:FloatByte}
    neg = hi <= 0
    if neg
        lo, hi = -hi, -lo                       # exact negation
    end
    (lo < 0 || !(hi > 0)) && return nothing     # straddles zero: refine
    # Overflow guard: everything at or above 2^(emax+2) saturates the same
    # way under every mode; synthesize an above-all-boundaries carrier.
    if !isfinite(hi) || lo >= ldexp(BigFloat(1), emax(Fr) + 2)
        isfinite(lo) || return nothing
        lo >= ldexp(BigFloat(1), emax(Fr) + 2) || return nothing
        return project_code(Fr, ρ, Aligned(UInt64(1), Int32(emax(Fr) + 3), neg, true), e)
    end
    # Tiny guard: everything strictly between 0 and quantum·2^-40 rounds
    # identically; synthesize a below-all-boundaries carrier.
    tiny_exp = (1 - bias(Fr)) - pbits(Fr) - 40
    if hi <= ldexp(BigFloat(1), tiny_exp)
        return project_code(Fr, ρ, Aligned(UInt64(1), Int32(tiny_exp - 2), neg, true), e)
    end
    ex = exponent(hi)
    s = 61 - ex
    L = trunc(BigInt, ldexp(lo, s))
    H = trunc(BigInt, ldexp(hi, s))
    L == H || return nothing
    return project_code(Fr, ρ, Aligned(UInt64(L), Int32(-s), neg, true), e)
end

function project_via_enclosure(::Type{Fr}, ρ::Projection, e::UInt64,
                               mkenc, opname::Symbol) where {Fr<:FloatByte}
    p = ELEM_P_START
    while p <= ELEM_P_CEILING
        lo, hi = mkenc(p)
        lo == hi && return _project_exact_big(Fr, ρ, e, lo)
        code = _try_window(Fr, ρ, e, lo, hi)
        code !== nothing && return code
        p *= 2
    end
    throw(InconclusiveProjectionError(opname, (), ELEM_P_CEILING))
end

@inline _datum_big(d::Datum) =
    ldexp(BigFloat(d.neg ? -Int64(d.sig) : Int64(d.sig)), Int(d.exp))

_direct_enclosure(f, d::Datum, p::Int) =
    setprecision(BigFloat, p) do
        x = _datum_big(d)
        lo = setrounding(() -> f(x), BigFloat, RoundDown)
        hi = setrounding(() -> f(x), BigFloat, RoundUp)
        (lo, hi)
    end

# f(x)/π for f(x) of known sign; interval quotient with directed endpoints.
_pi_div_enclosure(f, d::Datum, p::Int) =
    setprecision(BigFloat, p) do
        x = _datum_big(d)
        tlo = setrounding(() -> f(x), BigFloat, RoundDown)
        thi = setrounding(() -> f(x), BigFloat, RoundUp)
        (tlo == thi && iszero(tlo)) && return (tlo, thi)
        pilo = setrounding(() -> BigFloat(pi), BigFloat, RoundDown)
        pihi = setrounding(() -> BigFloat(pi), BigFloat, RoundUp)
        lo = setrounding(() -> tlo >= 0 ? tlo / pihi : tlo / pilo, BigFloat, RoundDown)
        hi = setrounding(() -> thi >= 0 ? thi / pilo : thi / pihi, BigFloat, RoundUp)
        (lo, hi)
    end

# ±π/2 and ±π constants through the same loop.
_half_pi_enclosure(neg::Bool, p::Int) =
    setprecision(BigFloat, p) do
        lo = setrounding(() -> BigFloat(pi) / 2, BigFloat, RoundDown)
        hi = setrounding(() -> BigFloat(pi) / 2, BigFloat, RoundUp)
        neg ? (-hi, -lo) : (lo, hi)
    end
_pi_enclosure(neg::Bool, p::Int) =
    setprecision(BigFloat, p) do
        lo = setrounding(() -> BigFloat(pi), BigFloat, RoundDown)
        hi = setrounding(() -> BigFloat(pi), BigFloat, RoundUp)
        neg ? (-hi, -lo) : (lo, hi)
    end

# Exact dyadic constants (draft rows whose results are dyadic).
_project_dyadic(::Type{Fr}, ρ, e, num::Int, ex::Int, neg::Bool) where {Fr<:FloatByte} =
    project_code(Fr, ρ, Aligned(UInt64(num), Int32(ex), neg, false), e)

# "Just below 1 in magnitude": for asymptotes whose gap to ±1 is far below
# every rounding boundary (Tanh at |x| ≥ 64, ExpMinusOne at x ≤ −64).
_project_near_one(::Type{Fr}, ρ, e, neg::Bool) where {Fr<:FloatByte} =
    project_code(Fr, ρ, Aligned((UInt64(1) << 62) - 1, Int32(-62), neg, true), e)

# Saturation threshold beyond which tanh/expm1 enclosures cannot converge.
const ASYMPTOTE_X = 64.0

# --- the unary transcendental kernels --------------------------------------
#
# Each entry: op => (mpfr call or custom, zero row, +∞ row, −∞ row, and any
# domain rows). Rows are implemented in the draft's order.

# result kinds for the special rows
# (:nan | :pinf | :ninf | :zero | :value)

function defined_code(::ExpOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_NINF && return 0x00
    d = k == CODE_ZERO ? Datum(UInt64(0), Int32(0), false) : unpack(F, c)
    d.sig == 0 && return project_dyadic_one(Fr, ρ, e)         # e^0 = 1
    return project_via_enclosure(Fr, ρ, e, p -> _direct_enclosure(exp, d, p), :Exp)
end
project_dyadic_one(::Type{Fr}, ρ, e) where {Fr<:FloatByte} =
    _project_dyadic(Fr, ρ, e, 1, 0, false)

function defined_code(::Exp2Op, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_NINF && return 0x00
    k == CODE_ZERO && return project_dyadic_one(Fr, ρ, e)
    return project_via_enclosure(Fr, ρ, e,
        p -> _direct_enclosure(exp2, unpack(F, c), p), :Exp2)
end

for (OP, FN, NAME) in ((:LogOp, :log, :Log), (:Log2Op, :log2, :Log2))
    @eval function defined_code(::$OP, ::Type{Fr}, ρ::Projection, e::UInt64,
                                x::F) where {Fr<:FloatByte, F<:FloatByte}
        c = codepoint(x); k = classify_code(F, c)
        k == CODE_NAN  && return nan_code(Fr)
        k == CODE_NINF && return nan_code(Fr)
        k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
        k == CODE_ZERO && return project_special(Fr, ρ, CODE_NINF)   # Log(0) → −∞
        d = unpack(F, c)
        d.neg && return nan_code(Fr)
        return project_via_enclosure(Fr, ρ, e,
            p -> _direct_enclosure($FN, d, p), $(QuoteNode(NAME)))
    end
end

function defined_code(::LogOnePlusOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_NINF && return nan_code(Fr)                # Log(1 + (−∞)) = Log(−∞)
    k == CODE_ZERO && return 0x00                        # log1p(0) = 0
    d = unpack(F, c)
    dv = decode(x)
    dv == -1 && return project_special(Fr, ρ, CODE_NINF) # Log(0) row
    dv < -1 && return nan_code(Fr)
    return project_via_enclosure(Fr, ρ, e, p -> _direct_enclosure(log1p, d, p), :LogOnePlus)
end

function defined_code(::ExpMinusOneOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_NINF && return _project_dyadic(Fr, ρ, e, 1, 0, true)   # 0 − 1 = −1
    k == CODE_ZERO && return 0x00
    d = unpack(F, c)
    dv = decode(x)
    dv <= -ASYMPTOTE_X && return _project_near_one(Fr, ρ, e, true)   # (−1, −1+ε)
    return project_via_enclosure(Fr, ρ, e, p -> _direct_enclosure(expm1, d, p), :ExpMinusOne)
end

# Sin, Cos, Tan, SinPi, CosPi: NaN at ±∞; finite through MPFR directly
# (MPFR sinpi/cospi are exact at half-integers, tanpi at quarter-integers).
for (OP, FN, NAME) in ((:SinOp, :sin, :Sin), (:CosOp, :cos, :Cos), (:TanOp, :tan, :Tan),
                       (:SinPiOp, :sinpi, :SinPi), (:CosPiOp, :cospi, :CosPi))
    @eval function defined_code(::$OP, ::Type{Fr}, ρ::Projection, e::UInt64,
                                x::F) where {Fr<:FloatByte, F<:FloatByte}
        c = codepoint(x); k = classify_code(F, c)
        k == CODE_NAN && return nan_code(Fr)
        (k == CODE_PINF || k == CODE_NINF) && return nan_code(Fr)
        k == CODE_ZERO &&
            return $(FN in (:cos, :cospi)) ? project_dyadic_one(Fr, ρ, e) : 0x00
        return project_via_enclosure(Fr, ρ, e,
            p -> _direct_enclosure($FN, unpack(F, c), p), $(QuoteNode(NAME)))
    end
end

function defined_code(::TanPiOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN && return nan_code(Fr)
    (k == CODE_PINF || k == CODE_NINF) && return nan_code(Fr)
    k == CODE_ZERO && return 0x00
    d = unpack(F, c)
    # Pole rows: cos(Xπ) = 0 at half-odd integers; the sign of sin(Xπ)
    # decides the infinity. X mod 2 == 1/2 → +∞; X mod 2 == 3/2 → −∞.
    if d.exp <= -1                                     # X has fractional part
        r = _mod2_times4(d)                            # 4·(|X| mod 2), or -1
        if r == 2                                       # |X| mod 2 == 1/2
            return project_special(Fr, ρ, d.neg ? CODE_NINF : CODE_PINF)
        elseif r == 6                                   # |X| mod 2 == 3/2
            return project_special(Fr, ρ, d.neg ? CODE_PINF : CODE_NINF)
        end
    elseif true
        # integer X: tan(Xπ) = 0 exactly — MPFR tanpi returns ±0.
    end
    return project_via_enclosure(Fr, ρ, e, p -> _direct_enclosure(tanpi, d, p), :TanPi)
end

# 4·(|X| mod 2) when it is an exact multiple of 1/4, else -1.
function _mod2_times4(d::Datum)
    ex = Int(d.exp)
    ex < -2 && begin
        # magnitude sig·2^ex: fractional resolution finer than 1/4 → check
        # whether the low bits below 1/4 vanish.
        sh = -2 - ex
        sh >= 64 && return d.sig == 0 ? 0 : -1
        (d.sig & ((UInt64(1) << sh) - 1)) != 0 && return -1
        return Int((d.sig >> sh) & 0x07)
    end
    return Int((d.sig << (ex + 2)) & 0x07)             # ex ∈ -2:-1 here or worse
end

for (OP, FN, NAME) in ((:SinhOp, :sinh, :Sinh), (:ArcSinhOp, :asinh, :ArcSinh))
    @eval function defined_code(::$OP, ::Type{Fr}, ρ::Projection, e::UInt64,
                                x::F) where {Fr<:FloatByte, F<:FloatByte}
        c = codepoint(x); k = classify_code(F, c)
        k == CODE_NAN  && return nan_code(Fr)
        k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
        k == CODE_NINF && return project_special(Fr, ρ, CODE_NINF)
        k == CODE_ZERO && return 0x00
        return project_via_enclosure(Fr, ρ, e,
            p -> _direct_enclosure($FN, unpack(F, c), p), $(QuoteNode(NAME)))
    end
end

function defined_code(::CoshOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN && return nan_code(Fr)
    (k == CODE_PINF || k == CODE_NINF) && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_ZERO && return project_dyadic_one(Fr, ρ, e)
    return project_via_enclosure(Fr, ρ, e,
        p -> _direct_enclosure(cosh, unpack(F, c), p), :Cosh)
end

function defined_code(::TanhOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_dyadic_one(Fr, ρ, e)            # exactly 1
    k == CODE_NINF && return _project_dyadic(Fr, ρ, e, 1, 0, true)   # exactly −1
    k == CODE_ZERO && return 0x00
    d = unpack(F, c)
    dv = decode(x)
    # Asymptote sidedness is symbolic: Tanh(large) is 1 − ε, never 1.
    abs(dv) >= ASYMPTOTE_X && return _project_near_one(Fr, ρ, e, d.neg)
    return project_via_enclosure(Fr, ρ, e, p -> _direct_enclosure(tanh, d, p), :Tanh)
end

function defined_code(::ArcCoshOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_NINF && return nan_code(Fr)
    k == CODE_ZERO && return nan_code(Fr)               # X < 1
    d = unpack(F, c)
    dv = decode(x)
    dv < 1 && return nan_code(Fr)
    dv == 1 && return 0x00                              # arccosh(1) = 0
    return project_via_enclosure(Fr, ρ, e, p -> _direct_enclosure(acosh, d, p), :ArcCosh)
end

function defined_code(::ArcTanhOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN && return nan_code(Fr)
    (k == CODE_PINF || k == CODE_NINF) && return nan_code(Fr)
    k == CODE_ZERO && return 0x00
    d = unpack(F, c)
    dv = decode(x)
    dv == 1  && return project_special(Fr, ρ, CODE_PINF)
    dv == -1 && return project_special(Fr, ρ, CODE_NINF)
    abs(dv) > 1 && return nan_code(Fr)
    return project_via_enclosure(Fr, ρ, e, p -> _direct_enclosure(atanh, d, p), :ArcTanh)
end

for (OP, FN, NAME) in ((:ArcSinOp, :asin, :ArcSin), (:ArcCosOp, :acos, :ArcCos))
    @eval function defined_code(::$OP, ::Type{Fr}, ρ::Projection, e::UInt64,
                                x::F) where {Fr<:FloatByte, F<:FloatByte}
        c = codepoint(x); k = classify_code(F, c)
        k == CODE_NAN && return nan_code(Fr)
        (k == CODE_PINF || k == CODE_NINF) && return nan_code(Fr)
        if k == CODE_ZERO
            return $(FN === :asin) ? 0x00 :
                   project_via_enclosure(Fr, ρ, e, p -> _half_pi_enclosure(false, p), :ArcCos)
        end
        d = unpack(F, c)
        abs(decode(x)) > 1 && return nan_code(Fr)
        return project_via_enclosure(Fr, ρ, e,
            p -> _direct_enclosure($FN, d, p), $(QuoteNode(NAME)))
    end
end

function defined_code(::ArcTanOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_via_enclosure(Fr, ρ, e, p -> _half_pi_enclosure(false, p), :ArcTan)
    k == CODE_NINF && return project_via_enclosure(Fr, ρ, e, p -> _half_pi_enclosure(true, p), :ArcTan)
    k == CODE_ZERO && return 0x00
    return project_via_enclosure(Fr, ρ, e,
        p -> _direct_enclosure(atan, unpack(F, c), p), :ArcTan)
end

function defined_code(::ArcSinPiOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN && return nan_code(Fr)
    (k == CODE_PINF || k == CODE_NINF) && return nan_code(Fr)
    k == CODE_ZERO && return 0x00
    d = unpack(F, c)
    dv = decode(x)
    abs(dv) > 1 && return nan_code(Fr)
    abs(dv) == 1 && return _project_dyadic(Fr, ρ, e, 1, -1, d.neg)   # ±1/2
    return project_via_enclosure(Fr, ρ, e, p -> _pi_div_enclosure(asin, d, p), :ArcSinPi)
end

function defined_code(::ArcCosPiOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN && return nan_code(Fr)
    (k == CODE_PINF || k == CODE_NINF) && return nan_code(Fr)
    k == CODE_ZERO && return _project_dyadic(Fr, ρ, e, 1, -1, false)  # 1/2
    d = unpack(F, c)
    dv = decode(x)
    abs(dv) > 1 && return nan_code(Fr)
    dv == 1  && return 0x00
    dv == -1 && return project_dyadic_one(Fr, ρ, e)
    return project_via_enclosure(Fr, ρ, e, p -> _pi_div_enclosure(acos, d, p), :ArcCosPi)
end

function defined_code(::ArcTanPiOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return _project_dyadic(Fr, ρ, e, 1, -1, false)  # 1/2
    k == CODE_NINF && return _project_dyadic(Fr, ρ, e, 1, -1, true)   # −1/2
    k == CODE_ZERO && return 0x00
    d = unpack(F, c)
    dv = decode(x)
    abs(dv) == 1 && return _project_dyadic(Fr, ρ, e, 1, -2, d.neg)    # ±1/4
    return project_via_enclosure(Fr, ρ, e, p -> _pi_div_enclosure(atan, d, p), :ArcTanPi)
end

function defined_code(::SoftplusOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x); k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_NINF && return 0x00
    softplus = t -> log1p(exp(t))       # monotone chain; directed per pass
    d = k == CODE_ZERO ? Datum(UInt64(0), Int32(0), false) : unpack(F, c)
    if d.sig == 0
        # softplus(0) = log 2, irrational; enclose it.
        return project_via_enclosure(Fr, ρ, e,
            p -> _direct_enclosure(t -> log(2 + 0 * t), Datum(UInt64(1), Int32(0), false), p),
            :Softplus)
    end
    return project_via_enclosure(Fr, ρ, e, p -> _direct_enclosure(softplus, d, p), :Softplus)
end

# --- ArcTan2, ArcTan2Pi (§4.10.15/.16) -------------------------------------

# Finite-quadrant enclosure: arctan(Y/X) (X > 0), ± π adjustment (X < 0),
# optionally divided by π. All steps are directed interval arithmetic.
function _atan2_enclosure(dy::Datum, dx::Datum, overpi::Bool, p::Int)
    setprecision(BigFloat, p) do
        yb = _datum_big(dy); xb = _datum_big(dx)
        qlo = setrounding(() -> yb / xb, BigFloat, RoundDown)
        qhi = setrounding(() -> yb / xb, BigFloat, RoundUp)
        tlo = setrounding(() -> atan(qlo), BigFloat, RoundDown)
        thi = setrounding(() -> atan(qhi), BigFloat, RoundUp)
        if dx.neg                              # X < 0: add ±π
            pilo = setrounding(() -> BigFloat(pi), BigFloat, RoundDown)
            pihi = setrounding(() -> BigFloat(pi), BigFloat, RoundUp)
            if dy.neg                          # Y < 0: arctan(Y/X) − π
                tlo, thi = setrounding(() -> tlo - pihi, BigFloat, RoundDown),
                           setrounding(() -> thi - pilo, BigFloat, RoundUp)
            else                               # Y > 0: arctan(Y/X) + π
                tlo, thi = setrounding(() -> tlo + pilo, BigFloat, RoundDown),
                           setrounding(() -> thi + pihi, BigFloat, RoundUp)
            end
        end
        overpi || return (tlo, thi)
        pilo = setrounding(() -> BigFloat(pi), BigFloat, RoundDown)
        pihi = setrounding(() -> BigFloat(pi), BigFloat, RoundUp)
        lo = setrounding(() -> tlo >= 0 ? tlo / pihi : tlo / pilo, BigFloat, RoundDown)
        hi = setrounding(() -> thi >= 0 ? thi / pilo : thi / pihi, BigFloat, RoundUp)
        (lo, hi)
    end
end

for (OP, overpi, NAME) in ((:ArcTan2Op, false, :ArcTan2), (:ArcTan2PiOp, true, :ArcTan2Pi))
    @eval function defined_code(::$OP, ::Type{Fr}, ρ::Projection, e::UInt64,
                                y::Fy, x::Fx) where {Fr<:FloatByte, Fy<:FloatByte, Fx<:FloatByte}
        cy = codepoint(y); cx = codepoint(x)
        ky = classify_code(Fy, cy); kx = classify_code(Fx, cx)
        (ky == CODE_NAN || kx == CODE_NAN) && return nan_code(Fr)
        yinf = ky == CODE_PINF || ky == CODE_NINF
        xinf = kx == CODE_PINF || kx == CODE_NINF
        (ky == CODE_ZERO && kx == CODE_ZERO) && return nan_code(Fr)
        (yinf && xinf) && return nan_code(Fr)
        # ∗, +∞ → 0 ; 0, X>0 → 0
        if kx == CODE_PINF || (ky == CODE_ZERO && kx == CODE_FINITE && !unpack(Fx, cx).neg)
            return 0x00
        end
        # +∞, ∗ → π/2 ; Y>0, 0 → π/2  (and the negative mirrors)
        yposhalf = ky == CODE_PINF ||
            (kx == CODE_ZERO && ky == CODE_FINITE && !unpack(Fy, cy).neg)
        yneghalf = ky == CODE_NINF ||
            (kx == CODE_ZERO && ky == CODE_FINITE && unpack(Fy, cy).neg)
        if yposhalf || yneghalf
            $(overpi) && return _project_dyadic(Fr, ρ, e, 1, -1, yneghalf)
            return project_via_enclosure(Fr, ρ, e,
                p -> _half_pi_enclosure(yneghalf, p), $(QuoteNode(NAME)))
        end
        # Y ≥ 0, −∞ → π ; Y < 0, −∞ → −π ; 0, X<0 → π
        if kx == CODE_NINF || (ky == CODE_ZERO && kx == CODE_FINITE && unpack(Fx, cx).neg)
            neg = ky == CODE_FINITE && kx == CODE_NINF ? unpack(Fy, cy).neg : false
            $(overpi) && return _project_dyadic(Fr, ρ, e, 1, 0, neg)
            return project_via_enclosure(Fr, ρ, e,
                p -> _pi_enclosure(neg, p), $(QuoteNode(NAME)))
        end
        dy = unpack(Fy, cy); dx = unpack(Fx, cx)
        return project_via_enclosure(Fr, ρ, e,
            p -> _atan2_enclosure(dy, dx, $(overpi), p), $(QuoteNode(NAME)))
    end
end
