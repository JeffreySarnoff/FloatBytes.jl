# operations/selection.jl — extrema family (§4.11), clamping, comparisons
# (§4.12), predicates and classification (§4.13), neighbours (§4.16).
# Comparisons and predicates stop after decode: every datum is an exact
# Float64 dyadic, so Float64 comparisons realize the draft rows exactly
# (single NaN, no −0). They produce no code and take no ρ or entropy.

# Selected-operand carrier: (class, format-erased value). Selection returns
# an operand (as an extended-real) and projects it once into Fr.
@inline function _project_operand(::Type{Fr}, ρ::Projection, e::UInt64,
                                  ::Type{F}, c::UInt8, k::UInt8) where {Fr<:FloatByte, F<:FloatByte}
    k == CODE_NAN  && return nan_code(Fr)
    k == CODE_PINF && return project_special(Fr, ρ, CODE_PINF)
    k == CODE_NINF && return project_special(Fr, ρ, CODE_NINF)
    k == CODE_ZERO && return 0x00
    return project_datum(Fr, ρ, unpack(F, c), e)
end

# The ten min/max variants share one core parameterized by three flags:
#   bymag  — compare |X| vs |Y| (ties broken by signed compare)
#   ismin  — pick the smaller
#   skipnan / skipinf — Number / Finite variants replace NaN (resp. skip
#   infinities) with the other operand per the draft tables.
for (OP, bymag, ismin, skipnan, skipinf) in (
        (:MinimumOp,                false, true,  false, false),
        (:MaximumOp,                false, false, false, false),
        (:MinimumNumberOp,          false, true,  true,  false),
        (:MaximumNumberOp,          false, false, true,  false),
        (:MinimumMagnitudeOp,       true,  true,  false, false),
        (:MaximumMagnitudeOp,       true,  false, false, false),
        (:MinimumMagnitudeNumberOp, true,  true,  true,  false),
        (:MaximumMagnitudeNumberOp, true,  false, true,  false),
        (:MinimumFiniteOp,          false, true,  true,  true),
        (:MaximumFiniteOp,          false, false, true,  true))
    @eval function defined_code(::$OP, ::Type{Fr}, ρ::Projection, e::UInt64,
                                x::Fx, y::Fy) where {Fr<:FloatByte, Fx<:FloatByte, Fy<:FloatByte}
        cx = codepoint(x); cy = codepoint(y)
        kx = classify_code(Fx, cx); ky = classify_code(Fy, cy)
        if $skipnan
            if kx == CODE_NAN
                return _project_operand(Fr, ρ, e, Fy, cy, ky)     # NaN,NaN → NaN too
            elseif ky == CODE_NAN
                return _project_operand(Fr, ρ, e, Fx, cx, kx)
            end
        else
            (kx == CODE_NAN || ky == CODE_NAN) && return nan_code(Fr)
        end
        xinf = kx == CODE_PINF || kx == CODE_NINF
        yinf = ky == CODE_PINF || ky == CODE_NINF
        if $skipinf
            if xinf && yinf
                if kx != ky
                    return _project_operand(Fr, ρ, e, Fx, cx,
                                            $ismin ? CODE_NINF : CODE_PINF)
                end
                return _project_operand(Fr, ρ, e, Fx, cx, kx)
            elseif xinf
                return _project_operand(Fr, ρ, e, Fy, cy, ky)
            elseif yinf
                return _project_operand(Fr, ρ, e, Fx, cx, kx)
            end
        elseif xinf || yinf
            if $bymag
                # MinimumMagnitude: equal infs pass through, mixed → −∞;
                # MaximumMagnitude: +∞ dominates, else −∞; single inf → other.
                if xinf && yinf
                    if kx == ky
                        return _project_operand(Fr, ρ, e, Fx, cx, kx)
                    end
                    k = $ismin ? CODE_NINF :
                        ((kx == CODE_PINF || ky == CODE_PINF) ? CODE_PINF : CODE_NINF)
                    return _project_operand(Fr, ρ, e, Fx, cx, k)
                end
                if $ismin
                    return xinf ? _project_operand(Fr, ρ, e, Fy, cy, ky) :
                                  _project_operand(Fr, ρ, e, Fx, cx, kx)
                else
                    return xinf ? _project_operand(Fr, ρ, e, Fx, cx, kx) :
                                  _project_operand(Fr, ρ, e, Fy, cy, ky)
                end
            else
                # Minimum: any −∞ → −∞; +∞ yields the other operand.
                if $ismin
                    (kx == CODE_NINF || ky == CODE_NINF) &&
                        return _project_operand(Fr, ρ, e, Fx, cx, CODE_NINF)
                    return kx == CODE_PINF ? _project_operand(Fr, ρ, e, Fy, cy, ky) :
                                             _project_operand(Fr, ρ, e, Fx, cx, kx)
                else
                    (kx == CODE_PINF || ky == CODE_PINF) &&
                        return _project_operand(Fr, ρ, e, Fx, cx, CODE_PINF)
                    return kx == CODE_NINF ? _project_operand(Fr, ρ, e, Fy, cy, ky) :
                                             _project_operand(Fr, ρ, e, Fx, cx, kx)
                end
            end
        end
        dx = decode(x); dy = decode(y)
        local takex::Bool
        if $bymag
            ax = abs(dx); ay = abs(dy)
            takex = ax != ay ? (($ismin) ? ax < ay : ax > ay) :
                    (($ismin) ? dx < dy : dx >= dy)
        else
            takex = $ismin ? dx < dy : !(dx < dy)
        end
        return takex ? _project_operand(Fr, ρ, e, Fx, cx, kx) :
                       _project_operand(Fr, ρ, e, Fy, cy, ky)
    end
end

# --- Clamp (§4.11.4), rows implemented literally in order ------------------

function defined_code(::ClampOp, ::Type{Fr}, ρ::Projection, e::UInt64,
                      x::Fx, lo::Fl, hi::Fh) where {Fr<:FloatByte, Fx<:FloatByte,
                                                    Fl<:FloatByte, Fh<:FloatByte}
    cx = codepoint(x); cl = codepoint(lo); ch = codepoint(hi)
    kx = classify_code(Fx, cx); kl = classify_code(Fl, cl); kh = classify_code(Fh, ch)
    (kx == CODE_NAN || kl == CODE_NAN || kh == CODE_NAN) && return nan_code(Fr)
    dl = decode(lo); dh = decode(hi); dx = decode(x)
    dl > dh && return nan_code(Fr)
    (kl == CODE_PINF && kh == CODE_PINF) && return _project_operand(Fr, ρ, e, Fh, ch, kh)
    (kl == CODE_NINF && kh == CODE_NINF) && return _project_operand(Fr, ρ, e, Fl, cl, kl)
    kh == CODE_NINF && return nan_code(Fr)
    kl == CODE_PINF && return nan_code(Fr)
    if kx == CODE_PINF
        kh == CODE_PINF && return _project_operand(Fr, ρ, e, Fx, cx, kx)
        return _project_operand(Fr, ρ, e, Fh, ch, kh)
    end
    if kx == CODE_NINF
        kl == CODE_NINF && return _project_operand(Fr, ρ, e, Fx, cx, kx)
        return _project_operand(Fr, ρ, e, Fl, cl, kl)
    end
    if kl == CODE_NINF
        (kh == CODE_PINF || dx <= dh) && return _project_operand(Fr, ρ, e, Fx, cx, kx)
        return _project_operand(Fr, ρ, e, Fh, ch, kh)
    end
    if kh == CODE_PINF
        dx <= dl && return _project_operand(Fr, ρ, e, Fl, cl, kl)
        return _project_operand(Fr, ρ, e, Fx, cx, kx)
    end
    dx <= dl && return _project_operand(Fr, ρ, e, Fl, cl, kl)
    dx >= dh && return _project_operand(Fr, ρ, e, Fh, ch, kh)
    return _project_operand(Fr, ρ, e, Fx, cx, kx)
end

# --- comparisons (§4.12) ---------------------------------------------------

compare_less(x::FloatByte, y::FloatByte)          = decode(x) < decode(y)
compare_less_equal(x::FloatByte, y::FloatByte)    = decode(x) <= decode(y)
compare_equal(x::FloatByte, y::FloatByte)         = decode(x) == decode(y)
compare_greater_equal(x::FloatByte, y::FloatByte) = decode(x) >= decode(y)
compare_greater(x::FloatByte, y::FloatByte)       = decode(x) > decode(y)

"""
    total_order(x, y) -> Bool

P3109 total order (§4.12.1): the single NaN compares as the most negative
value — NaN first, the opposite end from Float64 sorting.
"""
function total_order(x::FloatByte, y::FloatByte)
    is_nan_code(typeof(x), codepoint(x)) && return true
    is_nan_code(typeof(y), codepoint(y)) && return false
    return decode(x) <= decode(y)
end

# --- predicates and classification (§4.13) ---------------------------------

@enum FloatByteClass::UInt8 begin
    ClsNaN
    ClsNegativeInfinity
    ClsNegativeNormal
    ClsNegativeSubnormal
    ClsZero
    ClsPositiveSubnormal
    ClsPositiveNormal
    ClsPositiveInfinity
end

"""
    floatclass(x::FloatByte) -> FloatByteClass

The classifier operation `Class` (draft Table 2): every code point falls
into exactly one of the eight classes.
"""
function floatclass(x::F) where {F<:FloatByte}
    c = codepoint(x)
    k = classify_code(F, c)
    k == CODE_NAN  && return ClsNaN
    k == CODE_PINF && return ClsPositiveInfinity
    k == CODE_NINF && return ClsNegativeInfinity
    k == CODE_ZERO && return ClsZero
    S, E, neg = decode_fields(F, c)
    subnormal = S < (1 << (pbits(F) - 1))              # implies E == 1 - B
    if neg
        return subnormal ? ClsNegativeSubnormal : ClsNegativeNormal
    else
        return subnormal ? ClsPositiveSubnormal : ClsPositiveNormal
    end
end

is_zero(x::FloatByte)      = codepoint(x) == 0x00
is_one(x::F) where {F<:FloatByte} = codepoint(x) == one_code(F)
is_nan(x::F) where {F<:FloatByte} = is_nan_code(F, codepoint(x))
function is_infinite(x::F) where {F<:FloatByte}
    k = classify_code(F, codepoint(x))
    return k == CODE_PINF || k == CODE_NINF
end
is_finite(x::FloatByte) = !is_infinite(x) && !is_nan(x)
function is_sign_minus(x::F) where {F<:FloatByte}
    k = classify_code(F, codepoint(x))
    k == CODE_NINF && return true
    (k == CODE_NAN || k == CODE_PINF || k == CODE_ZERO) && return false
    return unpack(F, codepoint(x)).neg
end
function is_normal(x::F) where {F<:FloatByte}
    cls = floatclass(x)
    return cls == ClsNegativeNormal || cls == ClsPositiveNormal
end
function is_subnormal(x::F) where {F<:FloatByte}
    cls = floatclass(x)
    return cls == ClsNegativeSubnormal || cls == ClsPositiveSubnormal
end

# --- neighbours (§4.16): code-point walks with explicit terminal rows ------

"""
    nextgreaterthan(x) -> same format

`NextGreaterThan`: least value comparing greater, or NaN
(`NextGreaterThan(+Inf) = NaN`, unlike IEEE-754 `nextUp`).
"""
function nextgreaterthan(x::F) where {F<:FloatByte}
    c = codepoint(x)
    k = classify_code(F, c)
    k == CODE_NAN && return fromcode(F, nan_code(F))
    if isextended(F)
        k == CODE_PINF && return fromcode(F, nan_code(F))
        if issigned(F) && k == CODE_NINF
            return fromcode(F, maxfinite_code_signed(F, true))   # MinFinite
        end
    end
    c == maxfinite_code(F) &&
        return fromcode(F, isextended(F) ? pinf_code(F)::UInt8 : nan_code(F))
    if issigned(F) && c > (UInt8(1) << (kbits(F) - 1))           # negative finite
        c == (UInt8(1) << (kbits(F) - 1)) | 0x01 && return fromcode(F, 0x00)
        return fromcode(F, c - 0x01)
    end
    return fromcode(F, c + 0x01)
end

"""
    nextlessthan(x) -> same format

`NextLessThan`: greatest value comparing less, or NaN.
"""
function nextlessthan(x::F) where {F<:FloatByte}
    c = codepoint(x)
    k = classify_code(F, c)
    k == CODE_NAN && return fromcode(F, nan_code(F))
    if !issigned(F) && c == 0x00
        return fromcode(F, nan_code(F))
    end
    if issigned(F)
        minfin = maxfinite_code_signed(F, true)
        if isextended(F)
            k == CODE_NINF && return fromcode(F, nan_code(F))
            c == minfin && return fromcode(F, ninf_code(F)::UInt8)
        else
            c == minfin && return fromcode(F, nan_code(F))
        end
    end
    isextended(F) && k == CODE_PINF && return fromcode(F, maxfinite_code(F))
    if issigned(F)
        c == 0x00 && return fromcode(F, (UInt8(1) << (kbits(F) - 1)) | 0x01)
        if c > (UInt8(1) << (kbits(F) - 1))                      # negative finite
            return fromcode(F, c + 0x01)
        end
    end
    return fromcode(F, c - 0x01)
end
