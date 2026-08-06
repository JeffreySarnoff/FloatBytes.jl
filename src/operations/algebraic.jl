# operations/algebraic.jl — Add..Hypot kernels (§8). Kernel shape:
#   1. special rows — pure code-point tests in the draft's normative order
#   2. unpack survivors
#   3. one window routine from result_protocol.jl
#   4. one projection
# Never assume IEEE-754 host behavior: x/0 = NaN, Recip(0) = NaN,
# RSqrt(0) = NaN, one zero, one NaN.

const ONE_DATUM = Datum(UInt64(1), Int32(0), false)

# --- Convert (§4.9.1) ------------------------------------------------------

function defined_code(::ConvertOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x)
    k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_NINF && return project_special(Fr, ρ, CODE_NINF)
    k == CODE_ZERO && return 0x00
    return project_datum(Fr, ρ, unpack(F, c), e)
end

# --- Abs, Negate (§4.10.1), CopySign (§4.10.2) -----------------------------

function defined_code(::AbsOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x)
    k = classify_code(F, c)
    k == CODE_NAN && return nan_code(Fr)
    (k == CODE_PINF || k == CODE_NINF) && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_ZERO && return 0x00
    d = unpack(F, c)
    return project_datum(Fr, ρ, Datum(d.sig, d.exp, false), e)
end

function defined_code(::NegateOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x)
    k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_NINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_NINF)
    k == CODE_ZERO && return 0x00                      # -0 ≡ 0: one zero
    return project_datum(Fr, ρ, negate(unpack(F, c)), e)
end

function defined_code(::CopySignOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::Fx, y::Fy) where {Fr<:FloatByte, Fx<:FloatByte, Fy<:FloatByte}
    cx = codepoint(x); cy = codepoint(y)
    kx = classify_code(Fx, cx); ky = classify_code(Fy, cy)
    (kx == CODE_NAN || ky == CODE_NAN) && return nan_code(Fr)
    yneg = ky == CODE_NINF ? true :
           ky == CODE_PINF ? false :
           ky == CODE_ZERO ? false : unpack(Fy, cy).neg    # Y ≥ 0 → +
    if kx == CODE_PINF || kx == CODE_NINF
        return project_special(Fr, ρ, yneg ? CODE_NINF : CODE_PINF)
    end
    kx == CODE_ZERO && return 0x00
    d = unpack(Fx, cx)
    return project_datum(Fr, ρ, Datum(d.sig, d.exp, yneg), e)
end

# --- Add, Subtract (§4.10.3) -----------------------------------------------

@inline function _add_core(::Type{Fr}, ρ::Projection, e::UInt64,
                           ::Type{Fx}, cx::UInt8, kx::UInt8,
                           ::Type{Fy}, cy::UInt8, ky::UInt8,
                           negy::Bool) where {Fr<:FloatByte, Fx<:FloatByte, Fy<:FloatByte}
    (kx == CODE_NAN || ky == CODE_NAN) && return nan_code(Fr)
    # negy folds Subtract into Add by flipping y's sign, specials included.
    if negy
        ky = ky == CODE_PINF ? CODE_NINF : ky == CODE_NINF ? CODE_PINF : ky
    end
    xpinf = kx == CODE_PINF; xninf = kx == CODE_NINF
    ypinf = ky == CODE_PINF; yninf = ky == CODE_NINF
    ((xpinf && yninf) || (xninf && ypinf)) && return nan_code(Fr)
    (xpinf || ypinf) && return project_special(Fr, ρ, CODE_PINF)
    (xninf || yninf) && return project_special(Fr, ρ, CODE_NINF)
    if kx == CODE_ZERO
        ky == CODE_ZERO && return 0x00
        dy = unpack(Fy, cy)
        return project_datum(Fr, ρ, negy ? negate(dy) : dy, e)
    elseif ky == CODE_ZERO
        return project_datum(Fr, ρ, unpack(Fx, cx), e)
    end
    dx = unpack(Fx, cx)
    dy = unpack(Fy, cy)
    negy && (dy = negate(dy))
    return project_code(Fr, ρ, window_sum(dx, dy), e)
end

function defined_code(::AddOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::Fx, y::Fy) where {Fr<:FloatByte, Fx<:FloatByte, Fy<:FloatByte}
    cx = codepoint(x); cy = codepoint(y)
    return _add_core(Fr, ρ, e, Fx, cx, classify_code(Fx, cx),
                     Fy, cy, classify_code(Fy, cy), false)
end

function defined_code(::SubtractOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::Fx, y::Fy) where {Fr<:FloatByte, Fx<:FloatByte, Fy<:FloatByte}
    cx = codepoint(x); cy = codepoint(y)
    return _add_core(Fr, ρ, e, Fx, cx, classify_code(Fx, cx),
                     Fy, cy, classify_code(Fy, cy), true)
end

# --- Multiply (§4.10.4) ----------------------------------------------------

# Symbolic product class of two operand classes; also reports the result
# sign for infinite products. Used by Multiply and (via the documented
# Add∘Multiply equivalence) FMA.
@inline function _mul_symbolic(::Type{Fx}, cx::UInt8, kx::UInt8,
                               ::Type{Fy}, cy::UInt8, ky::UInt8) where {Fx, Fy}
    (kx == CODE_NAN || ky == CODE_NAN) && return (CODE_NAN, false)
    xinf = kx == CODE_PINF || kx == CODE_NINF
    yinf = ky == CODE_PINF || ky == CODE_NINF
    if xinf || yinf
        (kx == CODE_ZERO || ky == CODE_ZERO) && return (CODE_NAN, false)  # 0 × ∞
        xneg = kx == CODE_NINF ? true : kx == CODE_PINF ? false : unpack(Fx, cx).neg
        yneg = ky == CODE_NINF ? true : ky == CODE_PINF ? false : unpack(Fy, cy).neg
        return (xor(xneg, yneg) ? CODE_NINF : CODE_PINF, false)
    end
    (kx == CODE_ZERO || ky == CODE_ZERO) && return (CODE_ZERO, false)
    return (CODE_FINITE, false)
end

function defined_code(::MultiplyOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::Fx, y::Fy) where {Fr<:FloatByte, Fx<:FloatByte, Fy<:FloatByte}
    cx = codepoint(x); cy = codepoint(y)
    kx = classify_code(Fx, cx); ky = classify_code(Fy, cy)
    cls, _ = _mul_symbolic(Fx, cx, kx, Fy, cy, ky)
    cls == CODE_NAN  && return nan_code(Fr)
    cls == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    cls == CODE_NINF && return project_special(Fr, ρ, CODE_NINF)
    cls == CODE_ZERO && return 0x00
    return project_datum(Fr, ρ, window_mul(unpack(Fx, cx), unpack(Fy, cy)), e)
end

# --- Divide (§4.10.5), Recip (§4.10.8) -------------------------------------

function defined_code(::DivideOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::Fx, y::Fy) where {Fr<:FloatByte, Fx<:FloatByte, Fy<:FloatByte}
    cx = codepoint(x); cy = codepoint(y)
    kx = classify_code(Fx, cx); ky = classify_code(Fy, cy)
    (kx == CODE_NAN || ky == CODE_NAN) && return nan_code(Fr)
    xinf = kx == CODE_PINF || kx == CODE_NINF
    yinf = ky == CODE_PINF || ky == CODE_NINF
    (xinf && yinf) && return nan_code(Fr)              # ±∞ / ±∞
    ky == CODE_ZERO && return nan_code(Fr)             # x/0 = NaN (draft NOTE 2)
    if xinf
        yneg = unpack(Fy, cy).neg                      # y finite nonzero here
        xneg = kx == CODE_NINF
        return project_special(Fr, ρ, xor(xneg, yneg) ? CODE_NINF : CODE_PINF)
    end
    yinf && return 0x00                                # x/±∞ = 0 (draft NOTE 1)
    kx == CODE_ZERO && return 0x00
    return project_code(Fr, ρ, window_div(unpack(Fx, cx), unpack(Fy, cy)), e)
end

function defined_code(::RecipOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x)
    k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_ZERO && return nan_code(Fr)              # Recip(0) = NaN
    (k == CODE_PINF || k == CODE_NINF) && return 0x00  # Recip(±∞) = 0
    return project_code(Fr, ρ, window_div(ONE_DATUM, unpack(F, c)), e)
end

# --- FMA (§4.10.6): specials via the documented Add(Multiply(x,y), z) rows -

function defined_code(::FMAOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::Fx, y::Fy, z::Fz) where {Fr<:FloatByte, Fx<:FloatByte,
                                                  Fy<:FloatByte, Fz<:FloatByte}
    cx = codepoint(x); cy = codepoint(y); cz = codepoint(z)
    kx = classify_code(Fx, cx); ky = classify_code(Fy, cy); kz = classify_code(Fz, cz)
    kz == CODE_NAN && return nan_code(Fr)
    pcls, _ = _mul_symbolic(Fx, cx, kx, Fy, cy, ky)
    pcls == CODE_NAN && return nan_code(Fr)
    ppinf = pcls == CODE_PINF; pninf = pcls == CODE_NINF
    zpinf = kz == CODE_PINF;   zninf = kz == CODE_NINF
    ((ppinf && zninf) || (pninf && zpinf)) && return nan_code(Fr)
    (ppinf || zpinf) && return project_special(Fr, ρ, CODE_PINF)
    (pninf || zninf) && return project_special(Fr, ρ, CODE_NINF)
    if pcls == CODE_ZERO
        kz == CODE_ZERO && return 0x00
        return project_datum(Fr, ρ, unpack(Fz, cz), e)
    elseif kz == CODE_ZERO
        return project_datum(Fr, ρ, window_mul(unpack(Fx, cx), unpack(Fy, cy)), e)
    end
    return project_code(Fr, ρ,
        window_fma(unpack(Fx, cx), unpack(Fy, cy), unpack(Fz, cz)), e)
end

# --- FAA (§4.10.7): three terms, routed through LimbAcc (§4.4) -------------

# Fixed limb count covering every K ≤ 8 ternary signature: worst span is
# emax(Binary8p1uf) − (2 − 128 − 1) + sig/carry slack < 5·64 bits.
const FAA_LIMBS = 5

@inline _emin_frame(::Type{F}) where {F<:FloatByte} = 2 - bias(F) - pbits(F)

function defined_code(::FAAOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::Fx, y::Fy, z::Fz) where {Fr<:FloatByte, Fx<:FloatByte,
                                                  Fy<:FloatByte, Fz<:FloatByte}
    cx = codepoint(x); cy = codepoint(y); cz = codepoint(z)
    kx = classify_code(Fx, cx); ky = classify_code(Fy, cy); kz = classify_code(Fz, cz)
    (kx == CODE_NAN || ky == CODE_NAN || kz == CODE_NAN) && return nan_code(Fr)
    anyp = kx == CODE_PINF || ky == CODE_PINF || kz == CODE_PINF
    anyn = kx == CODE_NINF || ky == CODE_NINF || kz == CODE_NINF
    (anyp && anyn) && return nan_code(Fr)
    anyp && return project_special(Fr, ρ, CODE_PINF)
    anyn && return project_special(Fr, ρ, CODE_NINF)
    emin = min(_emin_frame(Fx), _emin_frame(Fy), _emin_frame(Fz))
    acc = LimbAcc{FAA_LIMBS}()
    kx == CODE_FINITE && (acc = update(acc, unpack(Fx, cx), emin))
    ky == CODE_FINITE && (acc = update(acc, unpack(Fy, cy), emin))
    kz == CODE_FINITE && (acc = update(acc, unpack(Fz, cz), emin))
    return project_code(Fr, ρ, close_to_aligned(acc, emin), e)
end

# --- Sqrt, RSqrt (§4.10.8) -------------------------------------------------

function defined_code(::SqrtOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x)
    k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_NINF && return nan_code(Fr)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_ZERO && return 0x00
    d = unpack(F, c)
    d.neg && return nan_code(Fr)                       # X < 0 → NaN
    return project_code(Fr, ρ, window_sqrt(d), e)
end

function defined_code(::RSqrtOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::F) where {Fr<:FloatByte, F<:FloatByte}
    c = codepoint(x)
    k = classify_code(F, c)
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_NINF && return nan_code(Fr)
    k == CODE_ZERO && return nan_code(Fr)              # X ≤ 0 → NaN
    k == CODE_PINF && return 0x00                      # RSqrt(+∞) = 0
    d = unpack(F, c)
    d.neg && return nan_code(Fr)
    return project_code(Fr, ρ, window_rsqrt(d), e)
end

# --- Hypot (§4.10.14) ------------------------------------------------------

function defined_code(::HypotOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::Fx, y::Fy) where {Fr<:FloatByte, Fx<:FloatByte, Fy<:FloatByte}
    cx = codepoint(x); cy = codepoint(y)
    kx = classify_code(Fx, cx); ky = classify_code(Fy, cy)
    # Draft row order: NaN rows precede the infinity rows here.
    (kx == CODE_NAN || ky == CODE_NAN) && return nan_code(Fr)
    (kx == CODE_PINF || kx == CODE_NINF || ky == CODE_PINF || ky == CODE_NINF) &&
        return project_special(Fr, ρ, CODE_PINF)
    if kx == CODE_ZERO
        ky == CODE_ZERO && return 0x00
        dy = unpack(Fy, cy)
        return project_datum(Fr, ρ, Datum(dy.sig, dy.exp, false), e)
    elseif ky == CODE_ZERO
        dx = unpack(Fx, cx)
        return project_datum(Fr, ρ, Datum(dx.sig, dx.exp, false), e)
    end
    return project_code(Fr, ρ, window_hypot(unpack(Fx, cx), unpack(Fy, cy)), e)
end
