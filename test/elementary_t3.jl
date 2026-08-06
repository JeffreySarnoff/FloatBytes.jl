# T3 (elementary) — interval-checked transcendentals (§17.1).
# Independent check: recompute a directed enclosure at high precision in
# the TEST, convert both dyadic endpoints exactly to rationals, project
# both through the reference rational projector; when the endpoint codes
# agree, the production result must equal them.

using FloatBytes: kbits, fromcode, codepoint, decode, defined_code

const ELEM_OPS = [
    (Exp, exp), (Exp2, exp2), (Log, log), (Log2, log2),
    (LogOnePlus, log1p), (ExpMinusOne, expm1),
    (Sin, sin), (Cos, cos), (Tan, tan),
    (ArcSin, asin), (ArcCos, acos), (ArcTan, atan),
    (Sinh, sinh), (Cosh, cosh), (Tanh, tanh),
    (ArcSinh, asinh), (ArcCosh, acosh), (ArcTanh, atanh),
    (SinPi, sinpi), (CosPi, cospi), (TanPi, tanpi),
    (Softplus, t -> log1p(exp(t))),
]
const PI_OPS = [(ArcSinPi, asin), (ArcCosPi, acos), (ArcTanPi, atan)]

# Directed enclosure of f(v) at precision p, as exact rationals; nothing if
# f is undefined or infinite at v under host evaluation.
function test_enclosure(f, v::Float64, p::Int; overpi::Bool = false)
    setprecision(BigFloat, p) do
        x = BigFloat(v)                       # exact
        local lo, hi
        try
            lo = setrounding(() -> f(x), BigFloat, RoundDown)
            hi = setrounding(() -> f(x), BigFloat, RoundUp)
        catch
            return nothing
        end
        (isnan(lo) || isnan(hi)) && return nothing
        (isinf(lo) || isinf(hi)) && return nothing
        if overpi
            pilo = setrounding(() -> BigFloat(pi), BigFloat, RoundDown)
            pihi = setrounding(() -> BigFloat(pi), BigFloat, RoundUp)
            l2 = setrounding(() -> lo >= 0 ? lo / pihi : lo / pilo, BigFloat, RoundDown)
            h2 = setrounding(() -> hi >= 0 ? hi / pilo : hi / pihi, BigFloat, RoundUp)
            lo, hi = l2, h2
        end
        return (Rational{BigInt}(lo), Rational{BigInt}(hi))
    end
end

function check_elem(F, f, OP, hostf, c, ρ, rm, rs; overpi = false)
    x = fromcode(F, c)
    v = decode(x)
    got = defined_code(OP, F, ρ, UInt64(0), x)
    (isnan(v) || isinf(v)) && return             # special rows pinned elsewhere
    enc = test_enclosure(hostf, v, 512; overpi = overpi)
    enc === nothing && return
    rlo, rhi = enc
    clo = DR.project_ref(f, rm, rs, rlo)
    chi = DR.project_ref(f, rm, rs, rhi)
    clo == chi || return                          # undecided at 512 bits: skip
    got == UInt8(clo) || @info "elem mismatch" F OP c v got clo
    @test got == UInt8(clo)
end

@testset "T3 elementary" begin
    rhos = [(Projection(NearestTiesToEven(), SatNone()), :nte, :satnone),
            (Projection(TowardPositive(), SatNone()), :tp, :satnone)]
    for F in (Binary8p3se, Binary4p2sf, Binary6p2ue)
        f = rformat(F)
        for (ρ, rm, rs) in rhos
            for (OP, hostf) in ELEM_OPS, c in 0:2^kbits(F)-1
                check_elem(F, f, OP, hostf, c, ρ, rm, rs)
            end
            for (OP, hostf) in PI_OPS, c in 0:2^kbits(F)-1
                check_elem(F, f, OP, hostf, c, ρ, rm, rs; overpi = true)
            end
        end
    end

    # exact points and asymptote rows
    F = Binary8p4se
    ρ = Projection(NearestTiesToEven(), SatNone())
    one_ = one(F); zero_ = zero(F)
    pinf = fromcode(F, FloatBytes.pinf_code(F))
    ninf = fromcode(F, FloatBytes.ninf_code(F))
    @test decode(exp(zero_)) == 1.0
    @test iszero(log(one_))
    @test decode(apply(Exp2, F, ρ, F(3.0))) == 8.0
    @test decode(apply(Log2, F, ρ, F(8.0))) == 3.0
    @test iszero(sin(zero_)) && decode(cos(zero_)) == 1.0
    @test decode(apply(Tanh, F, ρ, pinf)) == 1.0
    @test decode(apply(Tanh, F, ρ, ninf)) == -1.0
    # Asymptote sidedness: the VALUE is 1 − ε, never 1. NTE correctly
    # projects it to 1; directed modes toward zero must stay strictly below.
    @test decode(tanh(maxfinite(F))) == 1.0
    let ρz = Projection(TowardZero(), SatNone())
        @test decode(apply(Tanh, F, ρz, maxfinite(F))) ==
              decode(nextlessthan(one(F)))
        @test decode(apply(ExpMinusOne, F, ρz, minfinite(F))) ==
              decode(nextgreaterthan(F(-1.0)))
    end
    @test decode(apply(ExpMinusOne, F, ρ, ninf)) == -1.0
    @test decode(apply(ExpMinusOne, F, ρ, minfinite(F))) == -1.0
    @test iszero(apply(Exp, F, ρ, ninf))
    @test isnan(apply(Log, F, ρ, F(-1.0)))
    @test decode(apply(Log, F, ρ, zero_)) == -Inf
    @test decode(apply(SinPi, F, ρ, F(0.5))) == 1.0
    @test iszero(apply(SinPi, F, ρ, F(1.0)))
    @test decode(apply(CosPi, F, ρ, F(1.0))) == -1.0
    @test decode(apply(TanPi, F, ρ, F(0.25))) == 1.0
    @test decode(apply(TanPi, F, ρ, F(0.5))) == Inf   # pole, sin > 0
    @test decode(apply(TanPi, F, ρ, F(1.5))) == -Inf  # pole, sin < 0
    @test decode(apply(ArcTanPi, F, ρ, one_)) == 0.25
    @test decode(apply(ArcTanPi, F, ρ, pinf)) == 0.5
    @test decode(apply(ArcSinPi, F, ρ, one_)) == 0.5
    @test decode(apply(ArcCosPi, F, ρ, F(-1.0))) == 1.0
    @test iszero(apply(ArcCosPi, F, ρ, one_))
    @test isnan(apply(ArcSin, F, ρ, F(2.0)))
    @test decode(apply(ArcTanh, F, ρ, one_)) == Inf
    @test isnan(apply(ArcCosh, F, ρ, F(0.5)))
    @test iszero(apply(ArcCosh, F, ρ, one_))

    # ArcTan2 quadrant rows + sampled finite pairs
    @test isnan(apply(ArcTan2, F, ρ, zero_, zero_))
    @test isnan(apply(ArcTan2, F, ρ, pinf, pinf))
    @test iszero(apply(ArcTan2, F, ρ, zero_, one_))
    @test decode(apply(ArcTan2Pi, F, ρ, one_, zero_)) == 0.5
    @test decode(apply(ArcTan2Pi, F, ρ, zero_, F(-1.0))) == 1.0
    @test decode(apply(ArcTan2Pi, F, ρ, F(-1.0), ninf)) == -1.0
    let f8 = rformat(F)
        for _ in 1:200
            cy, cx = rand(0:255), rand(0:255)
            y, x = fromcode(F, cy), fromcode(F, cx)
            (isfinite(y) && isfinite(x) && !iszero(y) && !iszero(x)) || continue
            got = defined_code(ArcTan2, F, ρ, UInt64(0), y, x)
            enc = test_enclosure(t -> atan(t, BigFloat(decode(x))), decode(y), 512)
            enc === nothing && continue
            clo = DR.project_ref(f8, :nte, :satnone, enc[1])
            chi = DR.project_ref(f8, :nte, :satnone, enc[2])
            clo == chi && @test got == UInt8(clo)
        end
    end
end
