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

# Docstrings for the exported operation singletons; each names its draft
# clause. Use as `apply(Op, Fr, ρ, xs...)`.
for (nm, arity, clause, desc) in (
        (:Add, 2, "§4.10.3", "addition"),
        (:Subtract, 2, "§4.10.3", "subtraction"),
        (:Multiply, 2, "§4.10.4", "multiplication"),
        (:Divide, 2, "§4.10.5", "division (`x/0 = NaN`, `x/±Inf = 0`)"),
        (:FMA, 3, "§4.10.6", "fused multiply–add `x*y + z`, one projection"),
        (:FAA, 3, "§4.10.7", "fused add–add `x + y + z`, exactly associative, one projection"),
        (:Sqrt, 1, "§4.10.8", "square root (negative → NaN)"),
        (:Recip, 1, "§4.10.8", "reciprocal (`Recip(0) = NaN`, `Recip(±Inf) = 0`)"),
        (:RSqrt, 1, "§4.10.8", "reciprocal square root (`x ≤ 0` → NaN)"),
        (:Hypot, 2, "§4.10.14", "`√(x² + y²)` — algebraic, exact-window evaluated"),
        (:Abs, 1, "§4.10.1", "absolute value"),
        (:Negate, 1, "§4.10.1", "negation (in unsigned formats nonzero results follow the saturation rows)"),
        (:CopySign, 2, "§4.10.2", "magnitude of `x` with the sign of `y`"),
        (:Convert, 1, "§4.9.1", "format conversion (the only cross-format tier-1 operation)"),
        (:Minimum, 2, "§4.11.1", "minimum, NaN-propagating"),
        (:Maximum, 2, "§4.11.1", "maximum, NaN-propagating"),
        (:MinimumNumber, 2, "§4.11.1", "minimum; a NaN operand is replaced by the other operand"),
        (:MaximumNumber, 2, "§4.11.1", "maximum; a NaN operand is replaced by the other operand"),
        (:MinimumMagnitude, 2, "§4.11.2", "minimum by magnitude, signed tie-break"),
        (:MaximumMagnitude, 2, "§4.11.2", "maximum by magnitude, signed tie-break"),
        (:MinimumMagnitudeNumber, 2, "§4.11.2", "NaN-replacing minimum by magnitude"),
        (:MaximumMagnitudeNumber, 2, "§4.11.2", "NaN-replacing maximum by magnitude"),
        (:MinimumFinite, 2, "§4.11.3", "minimum preferring finite operands"),
        (:MaximumFinite, 2, "§4.11.3", "maximum preferring finite operands"),
        (:Clamp, 3, "§4.11.4", "clamp `x` to `[lo, hi]` (`lo > hi` → NaN)"),
        (:Exp, 1, "§4.10.9", "`e^x`"),
        (:Exp2, 1, "§4.10.9", "`2^x`"),
        (:Log, 1, "§4.10.9", "natural logarithm (`Log(0) = −∞`, negative → NaN)"),
        (:Log2, 1, "§4.10.9", "base-2 logarithm"),
        (:LogOnePlus, 1, "§4.10.9", "`log(1 + x)`"),
        (:ExpMinusOne, 1, "§4.10.9", "`e^x − 1`"),
        (:Sin, 1, "§4.10.10", "sine (±Inf → NaN)"),
        (:Cos, 1, "§4.10.10", "cosine (±Inf → NaN)"),
        (:Tan, 1, "§4.10.10", "tangent (±Inf → NaN)"),
        (:ArcSin, 1, "§4.10.10", "arcsine (|x| > 1 → NaN)"),
        (:ArcCos, 1, "§4.10.10", "arccosine (|x| > 1 → NaN)"),
        (:ArcTan, 1, "§4.10.10", "arctangent (`ArcTan(±Inf) = ±π/2`)"),
        (:Sinh, 1, "§4.10.11", "hyperbolic sine"),
        (:Cosh, 1, "§4.10.11", "hyperbolic cosine"),
        (:Tanh, 1, "§4.10.11", "hyperbolic tangent (`Tanh(±Inf) = ±1` exactly)"),
        (:ArcSinh, 1, "§4.10.11", "inverse hyperbolic sine"),
        (:ArcCosh, 1, "§4.10.11", "inverse hyperbolic cosine (x < 1 → NaN)"),
        (:ArcTanh, 1, "§4.10.11", "inverse hyperbolic tangent (`ArcTanh(±1) = ±∞`)"),
        (:SinPi, 1, "§4.10.12", "`sin(πx)`, exact at half-integers"),
        (:CosPi, 1, "§4.10.12", "`cos(πx)`, exact at half-integers"),
        (:TanPi, 1, "§4.10.12", "`tan(πx)`; poles at half-odd integers give `±∞`"),
        (:ArcSinPi, 1, "§4.10.12", "`arcsin(x)/π`"),
        (:ArcCosPi, 1, "§4.10.12", "`arccos(x)/π`"),
        (:ArcTanPi, 1, "§4.10.12", "`arctan(x)/π` (`±Inf → ±1/2`)"),
        (:Softplus, 1, "§4.10.13", "`log(1 + e^x)`"),
        (:ArcTan2, 2, "§4.10.15", "two-argument arctangent `(y, x)` (`ArcTan2(0,0) = NaN`)"),
        (:ArcTan2Pi, 2, "§4.10.16", "two-argument arctangent divided by π"))
    args = arity == 1 ? "x" : arity == 2 ? "x, y" : "x, y, z"
    @eval @doc $(
        "`" * string(nm) * "`: operation singleton — draft " * clause * ", " * desc *
        ". Apply with `apply(" * string(nm) * ", Fr, ρ, " * args * ")`."
    ) $nm
end

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
