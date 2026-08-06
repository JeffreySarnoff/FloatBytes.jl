# policies.jl — rounding/saturation singletons, ρ, entropy classes (§6).
# Policy is dispatch, never a branch.

abstract type RoundingMode end
"""`NearestTiesToEven()`: round to nearest; ties to the even code (draft §4.7.4; for `P == 1` formats evenness is defined on the code point). The package default."""
struct NearestTiesToEven <: RoundingMode end
"""`NearestTiesToAway()`: round to nearest; ties away from zero."""
struct NearestTiesToAway <: RoundingMode end
"""`TowardPositive()`: directed rounding toward `+∞`."""
struct TowardPositive    <: RoundingMode end
"""`TowardNegative()`: directed rounding toward `−∞`."""
struct TowardNegative    <: RoundingMode end
"""`TowardZero()`: directed rounding toward zero (truncation)."""
struct TowardZero        <: RoundingMode end
"""`ToOdd()`: round-to-odd — truncate, then force the kept code odd when inexact."""
struct ToOdd             <: RoundingMode end
"""`StochasticA(N)`: stochastic rounding, variant A of draft §4.7.4 — round away iff `⌊ν·2^N⌋ + R ≥ 2^N` for the supplied `N`-bit entropy value `R`. `N ∈ 1:63`; this build implements `N ≤ 32`."""
struct StochasticA{N}    <: RoundingMode end
"""`StochasticB(N)`: stochastic rounding, variant B — `⌊ν·2^(N+1)⌋ + 2R + 1 ≥ 2^(N+1)`."""
struct StochasticB{N}    <: RoundingMode end
"""`StochasticC(N)`: stochastic rounding, variant C — `RNITE(ν·2^N) + R ≥ 2^N`."""
struct StochasticC{N}    <: RoundingMode end

# Stateless singletons parameterized by N only; the random word travels
# separately (§12). N is validated at mode construction; N == 64 would need
# a two-word widemul and is not a defined mode.
for M in (:StochasticA, :StochasticB, :StochasticC)
    @eval function $M(n::Integer)
        1 <= n <= 63 || throw(ArgumentError("stochastic N must be in 1:63, got $n"))
        return $M{Int(n)}()
    end
end

abstract type SaturationMode end
"""`SatFinite()`: out-of-range results (including infinities) clamp to the extremal finite values."""
struct SatFinite    <: SaturationMode end
"""`SatPropagate()`: representable infinities propagate; finite overflow clamps to the extremal finite values."""
struct SatPropagate <: SaturationMode end
"""`SatNone()`: draft §4.7.5 default behavior — overflow goes to `±Inf` where representable (NaN in finite domains and for negative results in unsigned formats), except that a rounding direction pointing inward forces the extremal finite value."""
struct SatNone      <: SaturationMode end

"""
    Projection{R<:RoundingMode, S<:SaturationMode}

A projection specification ρ = (rounding mode, saturation mode).
"""
struct Projection{R<:RoundingMode, S<:SaturationMode} end
Projection(r::RoundingMode, s::SaturationMode) = Projection{typeof(r), typeof(s)}()

"""
    DEFAULT_PROJECTION

The immutable package default ρ: `(NearestTiesToEven, SatNone)`. Used by
the projecting constructors and the Base operator surface unless
overridden with [`with_projection`](@ref).
"""
const DEFAULT_PROJECTION = Projection{NearestTiesToEven, SatNone}()

# RoundOf / SatOf (draft §4.15)
round_of(::Projection{R,S}) where {R,S} = R()
sat_of(::Projection{R,S}) where {R,S} = S()

# Entropy class of a rounding mode — a trait, not a hierarchy.
abstract type EntropyClass end
struct Deterministic <: EntropyClass end
struct Stochastic    <: EntropyClass end

entropy_class(::Type{<:RoundingMode}) = Deterministic()
entropy_class(::Type{<:StochasticA})  = Stochastic()
entropy_class(::Type{<:StochasticB})  = Stochastic()
entropy_class(::Type{<:StochasticC})  = Stochastic()
entropy_class(::Type{Projection{R,S}}) where {R,S} = entropy_class(R)
entropy_class(ρ::Projection) = entropy_class(typeof(ρ))

is_stochastic(ρ::Projection) = entropy_class(ρ) isa Stochastic

_stochastic_bits(::Type{StochasticA{N}}) where {N} = N
_stochastic_bits(::Type{StochasticB{N}}) where {N} = N
_stochastic_bits(::Type{StochasticC{N}}) where {N} = N
_stochastic_bits(::Type{<:RoundingMode}) = 0
_stochastic_bits(::Type{Projection{R,S}}) where {R,S} = _stochastic_bits(R)

# --- error taxonomy (§20) --------------------------------------------------

"""
    InconclusiveProjectionError

The rigorous enclosure loop hit its precision ceiling undecided. The
destination is never mutated and no guess is written.
"""
struct InconclusiveProjectionError <: Exception
    op::Symbol
    operands::Tuple
    precision::Int
end

"""
    EntropyRequiredError

A stochastic ρ was requested without an entropy source; thrown at
preparation, never mid-loop.
"""
struct EntropyRequiredError <: Exception
    what::String
end

"""
    WorkspaceError

A caller-owned workspace is too small; thrown before any mutation of the
destination (required and provided sizes carried).
"""
struct WorkspaceError <: Exception
    required::Int
    provided::Int
end

"""
    AliasingError

Overlapping array views where the operation's semantics would differ.
"""
struct AliasingError <: Exception
    description::String
end
