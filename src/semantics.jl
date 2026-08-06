# semantics.jl — the seam (§10). `defined_code` is the internal seam:
# scalar wrappers, table builders, dense loops, and conformance all call it
# and nothing lower. `apply` is the public seam wrapping it.

abstract type Operation end

# Arithmetic (§8.2)
struct AddOp <: Operation end;       const Add = AddOp()
struct SubtractOp <: Operation end;  const Subtract = SubtractOp()
struct MultiplyOp <: Operation end;  const Multiply = MultiplyOp()
struct DivideOp <: Operation end;    const Divide = DivideOp()
struct FMAOp <: Operation end;       const FMA = FMAOp()
struct FAAOp <: Operation end;       const FAA = FAAOp()
struct SqrtOp <: Operation end;      const Sqrt = SqrtOp()
struct RecipOp <: Operation end;     const Recip = RecipOp()
struct RSqrtOp <: Operation end;     const RSqrt = RSqrtOp()
struct HypotOp <: Operation end;     const Hypot = HypotOp()
struct AbsOp <: Operation end;       const Abs = AbsOp()
struct NegateOp <: Operation end;    const Negate = NegateOp()
struct CopySignOp <: Operation end;  const CopySign = CopySignOp()
struct ConvertOp <: Operation end;   const Convert = ConvertOp()

# Selection (§8.3), clamp
struct MinimumOp <: Operation end;               const Minimum = MinimumOp()
struct MaximumOp <: Operation end;               const Maximum = MaximumOp()
struct MinimumNumberOp <: Operation end;         const MinimumNumber = MinimumNumberOp()
struct MaximumNumberOp <: Operation end;         const MaximumNumber = MaximumNumberOp()
struct MinimumMagnitudeOp <: Operation end;      const MinimumMagnitude = MinimumMagnitudeOp()
struct MaximumMagnitudeOp <: Operation end;      const MaximumMagnitude = MaximumMagnitudeOp()
struct MinimumMagnitudeNumberOp <: Operation end;const MinimumMagnitudeNumber = MinimumMagnitudeNumberOp()
struct MaximumMagnitudeNumberOp <: Operation end;const MaximumMagnitudeNumber = MaximumMagnitudeNumberOp()
struct MinimumFiniteOp <: Operation end;         const MinimumFinite = MinimumFiniteOp()
struct MaximumFiniteOp <: Operation end;         const MaximumFinite = MaximumFiniteOp()
struct ClampOp <: Operation end;                 const Clamp = ClampOp()

# Transcendentals (§9)
struct ExpOp <: Operation end;         const Exp = ExpOp()
struct Exp2Op <: Operation end;        const Exp2 = Exp2Op()
struct LogOp <: Operation end;         const Log = LogOp()
struct Log2Op <: Operation end;        const Log2 = Log2Op()
struct LogOnePlusOp <: Operation end;  const LogOnePlus = LogOnePlusOp()
struct ExpMinusOneOp <: Operation end; const ExpMinusOne = ExpMinusOneOp()
struct SinOp <: Operation end;         const Sin = SinOp()
struct CosOp <: Operation end;         const Cos = CosOp()
struct TanOp <: Operation end;         const Tan = TanOp()
struct ArcSinOp <: Operation end;      const ArcSin = ArcSinOp()
struct ArcCosOp <: Operation end;      const ArcCos = ArcCosOp()
struct ArcTanOp <: Operation end;      const ArcTan = ArcTanOp()
struct SinhOp <: Operation end;        const Sinh = SinhOp()
struct CoshOp <: Operation end;        const Cosh = CoshOp()
struct TanhOp <: Operation end;        const Tanh = TanhOp()
struct ArcSinhOp <: Operation end;     const ArcSinh = ArcSinhOp()
struct ArcCoshOp <: Operation end;     const ArcCosh = ArcCoshOp()
struct ArcTanhOp <: Operation end;     const ArcTanh = ArcTanhOp()
struct SinPiOp <: Operation end;       const SinPi = SinPiOp()
struct CosPiOp <: Operation end;       const CosPi = CosPiOp()
struct TanPiOp <: Operation end;       const TanPi = TanPiOp()
struct ArcSinPiOp <: Operation end;    const ArcSinPi = ArcSinPiOp()
struct ArcCosPiOp <: Operation end;    const ArcCosPi = ArcCosPiOp()
struct ArcTanPiOp <: Operation end;    const ArcTanPi = ArcTanPiOp()
struct SoftplusOp <: Operation end;    const Softplus = SoftplusOp()
struct ArcTan2Op <: Operation end;     const ArcTan2 = ArcTan2Op()
struct ArcTan2PiOp <: Operation end;   const ArcTan2Pi = ArcTan2PiOp()

"""
    defined_code(op, Fr, ρ, entropy::UInt64, xs::FloatByte...) -> UInt8

The internal seam: one semantic pipeline — code-point special rows in the
draft's order, exact auxiliary computation over the closed extended reals,
one final projection. A table entry *is* a call to this function.

Entropy contract: deterministic modes ignore the word; stochastic modes
consume exactly the word passed (its top `N` bits), including for exact
and special results. `entropy = 0` under a stochastic mode is legal
input, not a sentinel.
"""
function defined_code end

# --- the public seam -------------------------------------------------------

"""
    apply(op, Fr, ρ, xs...; entropy = nothing, index = 0) -> Fr

Public operation entry point. `Fr` is the result format, `ρ` the
projection specification. Operands are `FloatByte` values (heterogeneous
formats allowed). `entropy` is required for stochastic ρ: either a raw
`UInt64` word (top `N` bits are consumed) or an `IndexedEntropy` keyed by
`index`. An unsupported operand is a `MethodError`, never a silent
fall-through to host arithmetic.
"""
function apply(op::Operation, ::Type{Fr}, ρ::Projection, xs::FloatByte...;
               entropy = nothing, index::Integer = 0) where {Fr<:FloatByte}
    validate_format(Fr)
    word = _entropy_word_for(ρ, entropy, index)
    return fromcode(Fr, defined_code(op, Fr, ρ, word, xs...))
end

function _entropy_word_for(ρ::Projection, entropy, index)
    if !is_stochastic(ρ)
        return UInt64(0)
    end
    N = _stochastic_bits(typeof(ρ))
    N > 32 && throw(ArgumentError(
        "stochastic N > 32 is not implemented in this build (exact-residue " *
        "carriers of ImplementationPlan §4.5/§7.1 are a deferred stage)"))
    entropy === nothing && throw(EntropyRequiredError(
        "stochastic projection $(typeof(ρ)) requires an entropy word or IndexedEntropy"))
    entropy isa UInt64 && return entropy
    entropy isa IndexedEntropy &&
        return entropy_word(entropy, Int(index), SLOT_ELEMENT, Val(N)) << (64 - N)
    throw(ArgumentError("unsupported entropy source $(typeof(entropy))"))
end
