# T1 — test/lattice: all 13 296 codes vs the reference (§17.1).

using FloatBytes: kbits, fromcode, codepoint, isvalidcode, decode, decode32,
    FORMAT_ALIASES, nan_code, classify_code, CODE_NAN, CODE_PINF, CODE_NINF

# reference decode → Float64 (exact: every datum is a dyadic in range)
function refval(f, c)
    v = DR.decode_ref(f, c)
    v === :nan && return NaN
    v === :pinf && return Inf
    v === :ninf && return -Inf
    return Float64(v)          # exact for all covered formats (asserted below)
end

@testset "T1 lattice" begin
    for F in FORMAT_ALIASES
        f = rformat(F)
        K = kbits(F)
        for c in 0:(2^K - 1)
            x = fromcode(F, c)
            @test codepoint(x) == UInt8(c)
            rv = refval(f, c)
            pv = decode(x)
            @test isequal(pv, rv)
            @test isequal(Float64(decode32(x)), rv)   # Float32 atlas exactness
            # exactness of the Float64 representation itself
            v = DR.decode_ref(f, c)
            if v isa DR.Q
                @test DR.Q(pv) == v
            end
        end
        # range checking
        @test_throws ArgumentError fromcode(F, 2^K)
        @test_throws ArgumentError fromcode(F, -1)
        @test isvalidcode(F, 2^K - 1) && !isvalidcode(F, 2^K)
    end

    # encode ∘ decode round trip: exact re-projection reproduces the code
    # (finite codes; specials round-trip through their own paths)
    for F in FORMAT_ALIASES, c in 0:(2^kbits(F) - 1)
        x = fromcode(F, c)
        isfinite(x) || continue
        @test codepoint(quantize(F, decode(x))) == UInt8(c)
    end

    # class partition: every code exactly one class, matching decode
    for F in FORMAT_ALIASES, c in 0:(2^kbits(F) - 1)
        x = fromcode(F, c)
        cls = floatclass(x)
        v = decode(x)
        expected = isnan(v) ? ClsNaN :
                   v == Inf ? ClsPositiveInfinity :
                   v == -Inf ? ClsNegativeInfinity :
                   v == 0 ? ClsZero :
                   abs(v) >= decode(minnormal(F)) ?
                       (v < 0 ? ClsNegativeNormal : ClsPositiveNormal) :
                       (v < 0 ? ClsNegativeSubnormal : ClsPositiveSubnormal)
        @test cls == expected
    end

    # neighbour graph laws (§8.5): strict monotone steps, mutual inverses
    for F in FORMAT_ALIASES
        for c in 0:(2^kbits(F) - 1)
            x = fromcode(F, c)
            n = nextgreaterthan(x)
            if !isnan(x) && !isnan(n)
                @test decode(n) > decode(x)
                @test codepoint(nextlessthan(n)) == codepoint(x)
            end
            m = nextlessthan(x)
            if !isnan(x) && !isnan(m)
                @test decode(m) < decode(x)
                @test codepoint(nextgreaterthan(m)) == codepoint(x)
            end
        end
        # terminals
        isnan(nextgreaterthan(maxfinite(F))) ||
            @test decode(nextgreaterthan(maxfinite(F))) == Inf
        if FloatBytes.isextended(F)
            @test isnan(nextgreaterthan(fromcode(F, FloatBytes.pinf_code(F))))
        end
    end

    # total order: NaN first, strictly increasing decoded order
    for F in FORMAT_ALIASES
        codes = [fromcode(F, c) for c in 0:(2^kbits(F) - 1)]
        sorted = sort(codes)
        @test isnan(sorted[1])
        vals = [decode(x) for x in sorted[2:end]]
        @test issorted(vals) && allunique(vals)
        @test all(total_order(sorted[i], sorted[i + 1]) for i in 1:length(sorted) - 1)
    end

    # show never throws, including invalid bytes fabricated by reinterpret
    io = IOBuffer()
    for F in (Binary3p1sf, Binary5p3ue, Binary8p4se)
        for c in 0x00:0xff
            show(io, reinterpret(F, c))
        end
    end
    @test true
end
