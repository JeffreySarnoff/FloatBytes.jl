# T2 — projection vs the exact rational reference (§17.1).
# Values per format: every adjacent-datum gap edge, midpoint, midpoint ± a
# quarter gap, exact datums, beyond-max, below-min; every deterministic
# mode; the saturation cross-product; stochastic small-N full R grids.

using FloatBytes: kbits, fromcode, decode, FORMAT_ALIASES

# positive finite datum values of a format, ascending, as exact rationals
function posdatums(f)
    vals = DR.Q[]
    for c in DR.rcodes(f)
        v = DR.decode_ref(f, c)
        v isa DR.Q && v > 0 && push!(vals, v)
    end
    return sort!(vals)
end

function t2_values(f)
    vs = posdatums(f)
    out = DR.Q[]
    # gap grid (adjacent datums): edges, midpoint, midpoint ± quarter
    step = max(1, length(vs) ÷ 16)                  # subsample large formats
    for i in 1:step:length(vs) - 1
        a, b = vs[i], vs[i + 1]
        mid = (a + b) / 2
        q = (b - a) / 4
        push!(out, a, mid, mid - q, mid + q)
    end
    push!(out, vs[end])
    # saturation edges: beyond max, half-ulp above max, way beyond
    mx = vs[end]
    push!(out, mx + (mx - vs[end - 1]) / 2, 2 * mx, 1000 * mx)
    # below minimum positive: half of it, a quarter, way below
    mn = vs[1]
    push!(out, mn / 2, mn / 4, mn / 1000)
    return out
end

function check_projection(F, f, v::DR.Q, mode, rmode, sat, rsat; entropy = UInt64(0))
    ρ = Projection(mode, sat)
    got = FloatBytes.codepoint(quantize(F, v, ρ; entropy = entropy))
    want = DR.project_ref(f, rmode, rsat, v)
    if got != UInt8(want)
        @info "T2 mismatch" F v mode sat got want
    end
    @test got == UInt8(want)
end

@testset "T2 projection" begin
    formats = EXHAUSTIVE ? collect(FORMAT_ALIASES) :
        [F for F in FORMAT_ALIASES if kbits(F) <= 5]
    push!(formats, Binary8p4se, Binary8p1uf, Binary8p3uf, Binary8p5se)
    for F in formats
        f = rformat(F)
        for v0 in t2_values(f), v in (v0, -v0)
            for (m, rm) in DET_MODES, (s, rs) in SATS
                check_projection(F, f, v, m, rm, s, rs)
            end
        end
    end

    # stochastic: small N, ALL R, three variants, fractions around half
    sformats = [Binary4p2se, Binary4p3uf, Binary3p1sf, Binary5p2ue, Binary8p3se]
    for F in sformats
        f = rformat(F)
        vals = t2_values(f)
        for N in (1, 3), R in 0:(2^N - 1)
            word = UInt64(R) << (64 - N)
            for (mk, rk) in ((StochasticA(N), :stochA), (StochasticB(N), :stochB),
                             (StochasticC(N), :stochC))
                for v0 in vals[1:2:end], v in (v0, -v0)
                    ρ = Projection(mk, SatNone())
                    got = FloatBytes.codepoint(quantize(F, v, ρ; entropy = word))
                    want = DR.project_ref(f, (rk, N, R), :satnone, v)
                    @test got == UInt8(want)
                end
            end
        end
        # max covered N edge pins
        for R in (0, 2^32 - 1), v in (vals[1], vals[end ÷ 2])
            word = UInt64(R) << 32
            got = FloatBytes.codepoint(quantize(F, v, Projection(StochasticA(32), SatNone());
                                                entropy = word))
            want = DR.project_ref(f, (:stochA, 32, R), :satnone, v)
            @test got == UInt8(want)
        end
    end
end
