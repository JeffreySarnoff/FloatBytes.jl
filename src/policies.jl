# policies.jl — rounding/saturation singletons, ρ, entropy classes (§6).
# Policy is dispatch, never a branch.

abstract type RoundingMode end
struct NearestTiesToEven <: RoundingMode end
struct NearestTiesToAway <: RoundingMode end
struct TowardPositive    <: RoundingMode end
struct TowardNegative    <: RoundingMode end
struct TowardZero        <: RoundingMode end
struct ToOdd             <: RoundingMode end
struct StochasticA{N}    <: RoundingMode end     # N random bits, draft §4.7.4
struct StochasticB{N}    <: RoundingMode end
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
struct SatFinite    <: SaturationMode end
struct SatPropagate <: SaturationMode end
struct SatNone      <: SaturationMode end

"""
    Projection{R<:RoundingMode, S<:SaturationMode}

A projection specification ρ = (rounding mode, saturation mode).
"""
struct Projection{R<:RoundingMode, S<:SaturationMode} end
Projection(r::RoundingMode, s::SaturationMode) = Projection{typeof(r), typeof(s)}()

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

struct WorkspaceError <: Exception
    required::Int
    provided::Int
end

struct AliasingError <: Exception
    description::String
end
