# FloatBytes.jl — IEEE P3109 draft-1 byte-sized floating point for Julia.
#
# Implementation of docs/design/ImplementationPlan.md (stages 0–8 core;
# see docs/design/ImplementationStatus.md for stage coverage).
#
# Load order (single source of truth for the include sequence; a test
# greps that this comment matches the includes below):
#   result_protocol.jl        tier-1 window substrate; zero package deps
#   domain.jl                 FloatByte, traits, fromcode, decode atlas
#   formats.jl                the 120 generated aliases
#   show.jl                   display; loaded before anything that throws
#   policies.jl               rounding/saturation singletons, ρ, errors
#   projection.jl             round_to_precision ∘ saturate ∘ encode
#   execution/stochastic.jl   indexed entropy
#   semantics.jl              Operation singletons, defined_code, apply
#   operations/algebraic.jl   Add..Hypot window kernels
#   operations/selection.jl   min/max family, clamp, comparisons, class
#   operations/elementary.jl  transcendental enclosure loop (MPFR here only)
#   external.jl               Float64/Float32 bridges, quantize
#   julia_adapters.jl         deliberate Base surface
#
# Pinned draft revision: IEEE P3109 working draft D1 (docs/design/IEEE_D1.md).
# This package does not claim P3109 conformance; the draft is a draft.

module FloatBytes

import Random

include("result_protocol.jl")
include("domain.jl")
include("formats.jl")
include("show.jl")
include("policies.jl")
include("projection.jl")
include("execution/stochastic.jl")
include("semantics.jl")
include("operations/algebraic.jl")
include("operations/selection.jl")
include("operations/elementary.jl")
include("external.jl")
include("julia_adapters.jl")

const DRAFT_REVISION = "IEEE P3109 WG draft D1 (repository copy docs/design/IEEE_D1.md)"

# --- exports ---------------------------------------------------------------

export FloatByte, SignedFormat, UnsignedFormat, FiniteDomain, ExtendedDomain
# (the 120 aliases are exported from formats.jl)

export quantize, apply, with_projection, DEFAULT_PROJECTION
export Projection, NearestTiesToEven, NearestTiesToAway, TowardPositive,
       TowardNegative, TowardZero, ToOdd, StochasticA, StochasticB, StochasticC,
       SatFinite, SatPropagate, SatNone
export IndexedEntropy

# operation singletons
export Add, Subtract, Multiply, Divide, FMA, FAA, Sqrt, Recip, RSqrt, Hypot,
       Abs, Negate, CopySign, Convert,
       Minimum, Maximum, MinimumNumber, MaximumNumber,
       MinimumMagnitude, MaximumMagnitude,
       MinimumMagnitudeNumber, MaximumMagnitudeNumber,
       MinimumFinite, MaximumFinite, Clamp,
       Exp, Exp2, Log, Log2, LogOnePlus, ExpMinusOne,
       Sin, Cos, Tan, ArcSin, ArcCos, ArcTan,
       Sinh, Cosh, Tanh, ArcSinh, ArcCosh, ArcTanh,
       SinPi, CosPi, TanPi, ArcSinPi, ArcCosPi, ArcTanPi,
       Softplus, ArcTan2, ArcTan2Pi

# format-level and value queries
export maxfinite, minfinite, minpositive, maxsubnormal, minnormal,
       bitwidth, signedness, domain, exponent_bitwidth,
       trailing_significand_bits, exponent_bias

# predicates, classification, ordering, neighbours
export FloatByteClass, ClsNaN, ClsNegativeInfinity, ClsNegativeNormal,
       ClsNegativeSubnormal, ClsZero, ClsPositiveSubnormal,
       ClsPositiveNormal, ClsPositiveInfinity
export floatclass, total_order, nextgreaterthan, nextlessthan
export compare_less, compare_less_equal, compare_equal,
       compare_greater_equal, compare_greater

# errors
export InconclusiveProjectionError, EntropyRequiredError, WorkspaceError,
       AliasingError

# stable-but-unexported names (§3.2 API surface)
public fromcode, codepoint, isvalidcode, decode, format_id, format_params,
       format_type, FORMAT_ALIASES, entropy_word, defined_code, DRAFT_REVISION

end # module
