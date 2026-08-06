using Test
using FloatBytes

# The independent reference authority (§15): zero imports from src/.
include(joinpath(@__DIR__, "..", "reference", "draft_model.jl"))
using .DraftReference
const DR = DraftReference

# Tier control: FLOATBYTES_EXHAUSTIVE=1 widens the enumerations.
const EXHAUSTIVE = get(ENV, "FLOATBYTES_EXHAUSTIVE", "0") == "1"

# mode/sat mapping between production singletons and reference symbols
const DET_MODES = [
    (NearestTiesToEven(), :nte), (NearestTiesToAway(), :nta),
    (TowardPositive(), :tp), (TowardNegative(), :tn),
    (TowardZero(), :tz), (ToOdd(), :toodd),
]
const SATS = [(SatNone(), :satnone), (SatFinite(), :satfinite),
              (SatPropagate(), :satpropagate)]

rformat(F) = DR.RFormat(FloatBytes.kbits(F), FloatBytes.pbits(F),
                        FloatBytes.issigned(F), FloatBytes.isextended(F))

@testset "FloatBytes" begin
    include("stage0_envelope.jl")
    include("lattice.jl")
    include("projection_t2.jl")
    include("semantics_t3.jl")
    include("elementary_t3.jl")
    include("stochastic_t5.jl")
    include("julia_interface_t8.jl")
    include("structural.jl")
end
