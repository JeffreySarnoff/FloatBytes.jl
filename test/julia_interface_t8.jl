# T8 — Julia surface (§17.1): inference, generic algorithms, hashing,
# promotion, scoped projection, rand, parse, allocation of hot ops.

using FloatBytes: codepoint, fromcode, decode
using Random: Xoshiro

@testset "T8 julia interface" begin
    F = Binary8p4se
    x = F(1.5); y = F(2.25)

    @test @inferred(x + y) === F(3.75)
    @test @inferred(x * y) isa F
    @test @inferred(decode(x)) === 1.5
    @test @inferred(sqrt(y)) === F(1.5)
    @test @inferred(fma(x, y, y)) isa F

    # identities and queries
    @test iszero(zero(F)) && isone(one(F))
    @test decode(typemax(F)) == Inf && decode(typemin(F)) == -Inf
    @test decode(typemax(Binary8p4sf)) == decode(maxfinite(Binary8p4sf))
    @test iszero(typemin(Binary8p4uf))
    @test floatmax(F) === maxfinite(F)
    @test precision(F) == 4
    @test bitwidth(F) == 8
    @test exponent_bias(F) == 8
    @test isnan(maxsubnormal(Binary8p1uf))          # P == 1: no subnormals
    @test decode(eps(F)) == 2.0^-3

    # host conversions: exact both ways, one zero, one NaN
    @test Float64(x) == 1.5 && Float32(x) === 1.5f0
    @test BigFloat(x) == big"1.5"
    @test F(-0.0) === zero(F)                        # −0 collapses
    @test isnan(F(NaN)) && isnan(F(NaN32))
    @test decode(F(Inf)) == Inf
    @test F(Float32(1.5)) === x
    @test quantize(F, 1 // 3) isa F
    @test quantize(F, big"0.25") === F(0.25)
    @test F(3) === F(3.0)

    # promotion: Float64 target, exact
    @test x + 1.0 === 2.5
    @test 2.0f0 * x === 3.0f0
    @test_throws ArgumentError x + Binary8p3se(1.0) # cross-format: explicit apply

    # scoped projection override
    r1 = with_projection(() -> F(1.0) / F(3.0), Projection(TowardPositive(), SatNone()))
    r2 = with_projection(() -> F(1.0) / F(3.0), Projection(TowardNegative(), SatNone()))
    @test decode(r1) > 1 / 3 > decode(r2)

    # sorting: NaN first (P3109 total order), Dict/hash, Set
    rng = Xoshiro(42)
    v = rand(rng, F, 200)
    push!(v, fromcode(F, FloatBytes.nan_code(F)))
    s = sort(v)
    @test isnan(s[1])
    @test issorted(decode.(s[2:end]))
    d = Dict(zero(F) => 1, one(F) => 2)
    @test d[zero(F)] == 1 && d[one(F)] == 2
    @test length(Set([F(1.0), F(1.0), F(2.0)])) == 2

    # rand: finite datums only
    @test all(isfinite, rand(rng, F, 1000))
    @test all(isfinite, rand(rng, Binary3p1se, 200))

    # parse / show round trip on values
    @test parse(F, "1.5") === x
    @test tryparse(F, "zzz") === nothing
    @test occursin("Binary8p4se", sprint(show, x))
    @test occursin("1.5", sprint(show, x))

    # broadcast and zero-copy byte views
    a = [F(1.0), F(2.0), F(3.0)]
    @test decode.(a .+ F(1.0)) == [2.0, 3.0, 4.0]
    bytes = reinterpret(UInt8, a)
    @test bytes == codepoint.(a)
    @test sizeof(a) == 3 && isbitstype(F) && sizeof(F) == 1
    @test Base.elsize(Vector{F}) == 1

    # min/max/clamp draft semantics on the Base surface
    nan = fromcode(F, FloatBytes.nan_code(F))
    @test isnan(min(nan, x)) && isnan(max(x, nan))
    @test clamp(F(5.0), F(1.0), F(2.0)) === F(2.0)
    @test copysign(x, F(-1.0)) === -x
    @test abs(-x) === x

    # allocation: warm hot paths are allocation-free (W3 sample)
    f_add(a, b) = a + b
    f_add(x, y)
    @test @allocated(f_add(x, y)) == 0
    f_dec(a) = decode(a)
    f_dec(x)
    @test @allocated(f_dec(x)) == 0
end
