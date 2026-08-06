# julia_adapters.jl — the deliberate AbstractFloat surface (§14).
# Base operators require same-format operands under the scoped projection
# (default NTE, SatNone); cross-format arithmetic is a MethodError with
# `apply`/`Convert` as the documented path. Promotion targets Float64 —
# exact for every datum.

using Base.ScopedValues: ScopedValue, with

const PROJECTION_SCOPE = ScopedValue{Projection}(DEFAULT_PROJECTION)

"""
    with_projection(f, ρ)

Run `f()` with the Base-operator projection scope set to ρ. Outside the
performance interface: explicit `apply(op, Fr, ρ, ...)` is the fast path;
only the Base-operator adapters read this scope. No mutable global
numerical policy exists anywhere.
"""
with_projection(f, ρ::Projection) = with(f, PROJECTION_SCOPE => ρ)

# The scope read returns an abstractly-typed Projection; pinning the
# default gives Base operators a statically-dispatched, allocation-free
# fast path (the override path is documented as outside the performance
# interface, §14.2).
@inline function _scoped_call(op, ::Type{F}, xs...) where {F<:FloatByte}
    ρ = PROJECTION_SCOPE[]
    if ρ === DEFAULT_PROJECTION
        return fromcode(F, defined_code(op, F, DEFAULT_PROJECTION, UInt64(0), xs...))
    else
        return fromcode(F, defined_code(op, F, ρ, UInt64(0), xs...))
    end
end

# --- constructors and conversions ------------------------------------------

(::Type{F})(x::F) where {F<:FloatByte} = x
(::Type{Fr})(x::F) where {Fr<:FloatByte, F<:FloatByte} =
    _scoped_call(Convert, Fr, x)

Base.convert(::Type{F}, x::Real) where {F<:FloatByte} = F(x)

# --- identities, extrema, format-value queries (§14.1, draft §4.14) --------

Base.zero(::Type{F}) where {F<:FloatByte} = fromcode(F, 0x00)
Base.one(::Type{F}) where {F<:FloatByte} = fromcode(F, one_code(F))
Base.zero(x::FloatByte) = zero(typeof(x))
Base.one(x::FloatByte) = one(typeof(x))

"""maxfinite(F): `MaxFiniteOf` — maximum finite value of the format."""
maxfinite(::Type{F}) where {F<:FloatByte} = fromcode(F, maxfinite_code(F))
"""minfinite(F): `MinFiniteOf` — `-maxfinite` for signed formats, zero for unsigned."""
minfinite(::Type{F}) where {F<:FloatByte} =
    issigned(F) ? fromcode(F, maxfinite_code_signed(F, true)) : fromcode(F, 0x00)
"""minpositive(F): `MinPositiveOf` — least strictly positive value."""
minpositive(::Type{F}) where {F<:FloatByte} = fromcode(F, minpositive_code(F))
"""maxsubnormal(F): `MaxSubnormalOf`, or NaN when the format has no subnormals."""
maxsubnormal(::Type{F}) where {F<:FloatByte} =
    pbits(F) == 1 ? fromcode(F, nan_code(F)) :
                    fromcode(F, UInt8((1 << (pbits(F) - 1)) - 1))
"""minnormal(F): `MinNormalOf` — least positive normal value."""
minnormal(::Type{F}) where {F<:FloatByte} =
    fromcode(F, UInt8(1 << (pbits(F) - 1)))

Base.typemax(::Type{F}) where {F<:FloatByte} =
    isextended(F) ? fromcode(F, pinf_code(F)::UInt8) : maxfinite(F)
Base.typemin(::Type{F}) where {F<:FloatByte} =
    (issigned(F) && isextended(F)) ? fromcode(F, ninf_code(F)::UInt8) : minfinite(F)
Base.floatmax(::Type{F}) where {F<:FloatByte} = maxfinite(F)
Base.floatmin(::Type{F}) where {F<:FloatByte} = minnormal(F)

Base.precision(::Type{F}) where {F<:FloatByte} = pbits(F)
"""bitwidth(F): `BitwidthOf`."""
bitwidth(::Type{F}) where {F<:FloatByte} = kbits(F)
"""signedness(F): `SignednessOf` — `SignedFormat` or `UnsignedFormat`."""
signedness(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = S
"""domain(F): `DomainOf` — `FiniteDomain` or `ExtendedDomain`."""
domain(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = D
"""exponent_bitwidth(F): `ExponentBitwidthOf`."""
exponent_bitwidth(::Type{F}) where {F<:FloatByte} = expbits(F)
"""trailing_significand_bits(F): `TrailingSignificandBitwidthOf`."""
trailing_significand_bits(::Type{F}) where {F<:FloatByte} = trailingbits(F)
"""exponent_bias(F): `ExponentBiasOf`."""
exponent_bias(::Type{F}) where {F<:FloatByte} = bias(F)

function Base.eps(::Type{F}) where {F<:FloatByte}
    # ulp at one, per the Float64 convention.
    return quantize(F, exp2(1 - pbits(F)))
end
Base.eps(x::F) where {F<:FloatByte} = begin
    is_finite(x) || return fromcode(F, nan_code(F))
    n = nextgreaterthan(x)
    is_finite(n) ? F(abs(decode(n) - decode(x))) : eps(F)
end

# --- predicates ------------------------------------------------------------

Base.isnan(x::FloatByte) = is_nan(x)
Base.isinf(x::FloatByte) = is_infinite(x)
Base.isfinite(x::FloatByte) = is_finite(x)
Base.iszero(x::FloatByte) = is_zero(x)
Base.isone(x::FloatByte) = is_one(x)
Base.issubnormal(x::FloatByte) = is_subnormal(x)
Base.signbit(x::FloatByte) = is_sign_minus(x)

# --- comparison, ordering, hashing -----------------------------------------

Base.:(==)(x::FloatByte, y::FloatByte) = compare_equal(x, y)
Base.:(<)(x::FloatByte, y::FloatByte) = compare_less(x, y)
Base.:(<=)(x::FloatByte, y::FloatByte) = compare_less_equal(x, y)

# Julia total order for sorting: the P3109 total order — NaN sorts FIRST,
# inverting the Float64 convention (documented loudly).
function Base.isless(x::F, y::F) where {F<:FloatByte}
    xn = is_nan(x); yn = is_nan(y)
    xn && return !yn
    yn && return false
    return decode(x) < decode(y)
end

Base.hash(x::FloatByte, h::UInt) =
    hash((format_id(typeof(x)), codepoint(x)), hash(:FloatByte, h))

Base.cmp(x::F, y::F) where {F<:FloatByte} = isless(x, y) ? -1 : (isless(y, x) ? 1 : 0)

# --- promotion -------------------------------------------------------------

Base.promote_rule(::Type{F}, ::Type{Float64}) where {F<:FloatByte} = Float64
Base.promote_rule(::Type{F}, ::Type{Float32}) where {F<:FloatByte} = Float32
Base.promote_rule(::Type{F}, ::Type{Float16}) where {F<:FloatByte} = Float32
Base.promote_rule(::Type{F}, ::Type{<:Integer}) where {F<:FloatByte} = Float64
# No promote_rule between distinct FloatByte formats: a cross-format join
# is not pretended canonical; use apply/Convert.

# --- same-format Base operators under the scoped projection ----------------

for (op, OP) in ((:+, :Add), (:-, :Subtract), (:*, :Multiply), (:/, :Divide))
    @eval Base.$op(x::F, y::F) where {F<:FloatByte} = _scoped_call($OP, F, x, y)
    # Cross-format arithmetic is deliberately not defined: a cross-format
    # join is not pretended canonical. Point at the explicit seam.
    @eval Base.$op(x::FloatByte, y::FloatByte) = throw(ArgumentError(
        "cross-format $($(QuoteNode(op))) between $(typeof(x)) and $(typeof(y)) " *
        "is not defined; use apply($($(QuoteNode(OP))), Fr, ρ, x, y) or Convert"))
end
Base.:-(x::F) where {F<:FloatByte} =
    _scoped_call(Negate, F, x)
Base.abs(x::F) where {F<:FloatByte} =
    _scoped_call(Abs, F, x)
Base.sqrt(x::F) where {F<:FloatByte} =
    _scoped_call(Sqrt, F, x)
Base.inv(x::F) where {F<:FloatByte} =
    _scoped_call(Recip, F, x)
Base.fma(x::F, y::F, z::F) where {F<:FloatByte} =
    _scoped_call(FMA, F, x, y, z)
Base.muladd(x::F, y::F, z::F) where {F<:FloatByte} = fma(x, y, z)
Base.hypot(x::F, y::F) where {F<:FloatByte} =
    _scoped_call(Hypot, F, x, y)
Base.copysign(x::F, y::F) where {F<:FloatByte} =
    _scoped_call(CopySign, F, x, y)
Base.flipsign(x::F, y::F) where {F<:FloatByte} = signbit(y) ? -x : x

# Draft Minimum/Maximum semantics (NaN-propagating), not isless-based.
Base.min(x::F, y::F) where {F<:FloatByte} =
    _scoped_call(Minimum, F, x, y)
Base.max(x::F, y::F) where {F<:FloatByte} =
    _scoped_call(Maximum, F, x, y)
Base.clamp(x::F, lo::F, hi::F) where {F<:FloatByte} =
    _scoped_call(Clamp, F, x, lo, hi)

for (fn, OP) in ((:exp, :Exp), (:exp2, :Exp2), (:expm1, :ExpMinusOne),
                 (:log, :Log), (:log2, :Log2), (:log1p, :LogOnePlus),
                 (:sin, :Sin), (:cos, :Cos), (:tan, :Tan),
                 (:asin, :ArcSin), (:acos, :ArcCos), (:atan, :ArcTan),
                 (:sinh, :Sinh), (:cosh, :Cosh), (:tanh, :Tanh),
                 (:asinh, :ArcSinh), (:acosh, :ArcCosh), (:atanh, :ArcTanh),
                 (:sinpi, :SinPi), (:cospi, :CosPi))
    @eval Base.$fn(x::F) where {F<:FloatByte} = _scoped_call($OP, F, x)
end
Base.atan(y::F, x::F) where {F<:FloatByte} =
    _scoped_call(ArcTan2, F, y, x)

# --- parsing and rand ------------------------------------------------------

function Base.tryparse(::Type{F}, s::AbstractString) where {F<:FloatByte}
    b = tryparse(BigFloat, s)
    b === nothing && return nothing
    return quantize(F, b)
end
function Base.parse(::Type{F}, s::AbstractString) where {F<:FloatByte}
    r = tryparse(F, s)
    r === nothing && throw(ArgumentError("cannot parse $(repr(s)) as $(F)"))
    return r
end

# rand: uniform over finite datums (documented choice). The Sampler hook
# is required because FloatByte <: AbstractFloat would otherwise fall into
# Random's CloseOpen01 float machinery.
Random.Sampler(::Type{<:Random.AbstractRNG}, ::Type{F},
               ::Union{Val{1}, Val{Inf}}) where {F<:FloatByte} =
    Random.SamplerType{F}()
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{F}) where {F<:FloatByte}
    while true
        c = rand(rng, UInt8) & UInt8((1 << kbits(F)) - 1)
        k = classify_code(F, c)
        (k == CODE_NAN || k == CODE_PINF || k == CODE_NINF) && continue
        return fromcode(F, c)
    end
end
