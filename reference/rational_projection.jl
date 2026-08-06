# reference/rational_projection.jl — projection by exact rational
# quotient/remainder, transcribed literally from draft §4.7.3–.6.
# Decides every rounding question from exact arithmetic; never shares the
# production window/guard-bit code. Included at the end of draft_model.jl,
# inside module DraftReference.

# ρ: (mode, sat); mode ∈ :nte,:nta,:tp,:tn,:tz,:toodd or (:stochA|:stochB|:stochC, N, R)

# --- ωRoundToPrecision (draft §4.7.4) --------------------------------------
# Returns :zero | (S::BigInt, Qexp::Int, neg::Bool) | passthrough symbol.

function round_to_precision_ref(P::Int, B::Int, mode, X)
    X isa Symbol && return X                        # NaN, ±∞ pass through
    m = X.mag
    _is_zero(m) && return :zero
    e = floor_log2(m)
    Qexp = max(e, 1 - B) - P + 1
    Sfl = floor_scaled(m, -Qexp)                    # ⌊S̃⌋
    inexact = !is_int_scaled(m, -Qexp)
    away = round_away_ref(mode, m, Qexp, Sfl, inexact, X.neg, P, B)
    S = Sfl + (away ? 1 : 0)
    S == 0 && return :zero
    return (S, Qexp, X.neg)
end

function code_is_even_ref(Sfl::BigInt, Qexp::Int, P::Int, B::Int)
    P > 1 ? iseven(Sfl) : (Sfl == 0 || iseven(Qexp + B))
end

function round_away_ref(mode, m, Qexp, Sfl, inexact, neg, P, B)
    if mode === :tz
        return false
    elseif mode === :tp
        return inexact && !neg
    elseif mode === :tn
        return inexact && neg
    elseif mode === :nta
        # ν ≥ 1/2 ⟺ ⌊m·2^(1-Q)⌋ ≥ 2⌊S̃⌋ + 1
        return floor_scaled(m, 1 - Qexp) - 2 * Sfl >= 1
    elseif mode === :nte
        t1 = floor_scaled(m, 1 - Qexp) - 2 * Sfl
        nu_gt_half = t1 >= 1 && !is_int_scaled(m, 1 - Qexp)
        nu_eq_half = t1 >= 1 && is_int_scaled(m, 1 - Qexp)
        return nu_gt_half || (nu_eq_half && !code_is_even_ref(Sfl, Qexp, P, B))
    elseif mode === :toodd
        return inexact && code_is_even_ref(Sfl, Qexp, P, B)
    end
    kind, N, R = mode
    if kind === :stochA
        t = floor_scaled(m, N - Qexp) - Sfl * (BigInt(1) << N)
        return t + R >= BigInt(1) << N
    elseif kind === :stochB
        t = floor_scaled(m, N + 1 - Qexp) - Sfl * (BigInt(1) << (N + 1))
        return t + (2 * BigInt(R) + 1) >= BigInt(1) << (N + 1)
    elseif kind === :stochC
        M = BigInt(1) << N
        t = floor_scaled(m, N - Qexp) - Sfl * M
        t2 = floor_scaled(m, N + 1 - Qexp) - Sfl * 2 * M
        b = t2 - 2t                                 # next fraction bit
        rest = !is_int_scaled(m, N + 1 - Qexp)
        rn = t + ((b == 1 && (rest || isodd(t))) ? 1 : 0)
        return rn + R >= M
    end
    error("unknown mode $mode")
end

# --- ωSaturate (draft §4.7.5, rows in order) + ωEncode (§4.7.6) ------------

function encode_ref(f::RFormat, S::BigInt, Qexp::Int, neg::Bool)
    S == 0 && return 0
    P = f.P; B = rbias(f)
    # value = S·2^Q; re-derive draft-encode fields E, S', T
    v = Q(S) * Q(2)^Qexp
    E = max(floor_log2(v), 1 - B)
    Sd = floor_scaled(v, -(E + 1 - P))              # S' = v·2^{-E}·2^{P-1}
    T = Sd % (BigInt(1) << (P - 1))
    r = Sd < (BigInt(1) << (P - 1)) ? T : T + (BigInt(E + B) << (P - 1))
    return Int(r) + (neg ? 1 << (f.K - 1) : 0)
end

function project_ref(f::RFormat, mode, sat, X)
    X isa Q && (X = RVal(X))
    X === :nan && return rnan(f)
    R = round_to_precision_ref(f.P, rbias(f), mode, X)
    return saturate_encode_ref(f, mode, sat, R)
end

_dirmode(mode) = mode isa Symbol ? mode : :stoch

function saturate_encode_ref(f::RFormat, mode, sat, R)
    R === :zero && return 0
    mhi = decode_ref(f, rmaxfinite_code(f))::Q
    mlo = f.signed ? -mhi : Q(0)
    if R === :pinf || R === :ninf
        pos = R === :pinf
        if sat === :satfinite
            return pos ? rmaxfinite_code(f) :
                   (f.signed ? rmaxfinite_code(f) + (1 << (f.K - 1)) : 0)
        elseif sat === :satpropagate
            pos && return f.extended ? rpinf(f) : rmaxfinite_code(f)
            (f.signed && f.extended) && return rninf(f)
            return f.signed ? rmaxfinite_code(f) + (1 << (f.K - 1)) : 0
        else
            pos && return f.extended ? rpinf(f) : rnan(f)
            (f.signed && f.extended) && return rninf(f)
            return rnan(f)
        end
    end
    S, Qexp, neg = R
    v = Q(S) * Q(2)^Qexp * (neg ? -1 : 1)
    if mlo <= v <= mhi
        return encode_ref(f, S, Qexp, neg)
    end
    d = _dirmode(mode)
    if sat === :satfinite || sat === :satpropagate
        v > mhi && return rmaxfinite_code(f)
        return f.signed ? rmaxfinite_code(f) + (1 << (f.K - 1)) : 0
    end
    # SatNone
    if v > mhi && (d === :tz || d === :tn)
        return rmaxfinite_code(f)
    end
    if v < mlo && (d === :tz || d === :tp)
        return f.signed ? rmaxfinite_code(f) + (1 << (f.K - 1)) : 0
    end
    if f.extended
        v > mhi && return rpinf(f)
        f.signed && return rninf(f)
        return rnan(f)                              # negative in unsigned extended
    end
    return rnan(f)
end

# --- ω-operations on decoded values (draft §4.10 rows) ---------------------
# Values: :nan | :pinf | :ninf | Rational{BigInt} (signed). Sqrt-family
# finite results are RVal with SQ magnitude.

_isinf(v) = v === :pinf || v === :ninf
_val(v::Q) = RVal(v)

function omega_add(x, y)
    (x === :nan || y === :nan) && return :nan
    if _isinf(x) || _isinf(y)
        (x === :pinf && y === :ninf) && return :nan
        (x === :ninf && y === :pinf) && return :nan
        (x === :pinf || y === :pinf) && return :pinf
        return :ninf
    end
    return _val(x + y)
end

omega_subtract(x, y) =
    omega_add(x, y === :nan ? :nan : y === :pinf ? :ninf : y === :ninf ? :pinf : -y)

function omega_multiply(x, y)
    (x === :nan || y === :nan) && return :nan
    if _isinf(x) || _isinf(y)
        ((x isa Q && iszero(x)) || (y isa Q && iszero(y))) && return :nan
        sx = x === :pinf ? false : x === :ninf ? true : x < 0
        sy = y === :pinf ? false : y === :ninf ? true : y < 0
        return xor(sx, sy) ? :ninf : :pinf
    end
    return _val(x * y)
end

function omega_divide(x, y)
    (x === :nan || y === :nan) && return :nan
    (_isinf(x) && _isinf(y)) && return :nan
    (y isa Q && iszero(y)) && return :nan
    if _isinf(x)
        sx = x === :ninf
        sy = y < 0
        return xor(sx, sy) ? :ninf : :pinf
    end
    _isinf(y) && return _val(Q(0))
    return _val(x / y)
end

function omega_fma(x, y, z)
    (x === :nan || y === :nan || z === :nan) && return :nan
    p = omega_multiply(x, y)
    p === :nan && return :nan
    if _isinf(p) || _isinf(z)
        return omega_add(p, z)
    end
    return _val(x * y + z)                          # exact, single projection later
end

function omega_faa(x, y, z)
    (x === :nan || y === :nan || z === :nan) && return :nan
    anyp = x === :pinf || y === :pinf || z === :pinf
    anyn = x === :ninf || y === :ninf || z === :ninf
    (anyp && anyn) && return :nan
    anyp && return :pinf
    anyn && return :ninf
    return _val(x + y + z)
end

function omega_sqrt(x)
    x === :nan && return :nan
    x === :ninf && return :nan
    x === :pinf && return :pinf
    x < 0 && return :nan
    iszero(x) && return _val(Q(0))
    return RVal(false, SQ(x))
end

function omega_recip(x)
    x === :nan && return :nan
    _isinf(x) && return _val(Q(0))
    iszero(x) && return :nan
    return _val(1 / x)
end

function omega_rsqrt(x)
    x === :nan && return :nan
    x === :ninf && return :nan
    x === :pinf && return _val(Q(0))
    (x isa Q && x <= 0) && return :nan
    return RVal(false, SQ(1 / x))
end

function omega_hypot(x, y)
    (x === :nan || y === :nan) && return :nan       # NaN rows precede Inf rows
    (_isinf(x) || _isinf(y)) && return :pinf
    s = x * x + y * y
    iszero(s) && return _val(Q(0))
    return RVal(false, SQ(s))
end

omega_convert(x) = x isa Q ? _val(x) : x
