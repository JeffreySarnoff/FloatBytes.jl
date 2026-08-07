using Aqua
using JET

@testset "package quality" begin
    # FloatBytes deliberately extends Base's numerical interface, so method
    # ambiguity checking is covered by the focused interface tests instead.
    Aqua.test_all(FloatBytes; ambiguities = false)
    JET.test_package(FloatBytes; target_modules = (FloatBytes,))
end
