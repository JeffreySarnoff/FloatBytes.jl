# Stage 0 envelope (§19 stage 0): machine-check every §1 number.

using FloatBytes: kbits, pbits, issigned, isextended, bias, expbits, trailingbits,
    nan_code, pinf_code, ninf_code, one_code, maxfinite_code, format_id,
    format_params, format_type, FORMAT_ALIASES, FORMAT_OFFSETS, N_CODEPOINTS,
    fromcode, decode, DRAFT_REVISION

@testset "stage0 envelope" begin
    # 120 formats, 13 296 code points, 4K − 2 formats per K
    @test length(FORMAT_ALIASES) == 120
    @test N_CODEPOINTS == 13_296
    @test length(DR.all_formats()) == 120
    for K in 3:8
        n = count(F -> kbits(F) == K, collect(FORMAT_ALIASES))
        @test n == 4K - 2
    end

    # lattice id bijection + alias agreement + offsets
    for id in 1:120
        F = format_type(id)
        @test format_id(F) == id
        @test FORMAT_ALIASES[id] === F
        K, P, s, e = format_params(id)
        @test (kbits(F), pbits(F), issigned(F), isextended(F)) == (K, P, s, e)
        @test FORMAT_OFFSETS[id + 1] - FORMAT_OFFSETS[id] == 1 << K
    end

    # §1.1 table: bias, field widths, legal P ranges
    for F in FORMAT_ALIASES
        K, P = kbits(F), pbits(F)
        if issigned(F)
            @test 1 <= P <= K - 1
            @test bias(F) == 2^(K - P - 1)
            @test expbits(F) == K - P
        else
            @test 1 <= P <= K
            @test bias(F) == 2^(K - P)
            @test expbits(F) == K - P + 1
        end
        @test trailingbits(F) == P - 1
    end
    # maximum bias is 128 (K = 8, P = 1, unsigned)
    @test maximum(bias(F) for F in FORMAT_ALIASES) == 128
    @test bias(Binary8p1uf) == 128

    # §1.2 special code points
    for F in FORMAT_ALIASES
        K = kbits(F)
        @test nan_code(F) == UInt8(issigned(F) ? 2^(K - 1) : 2^K - 1)
        @test one_code(F) == UInt8(issigned(F) ? 2^(K - 2) : 2^(K - 1))
        if isextended(F)
            @test pinf_code(F) == UInt8(issigned(F) ? 2^(K - 1) - 1 : 2^K - 2)
        else
            @test pinf_code(F) === nothing
        end
        if isextended(F) && issigned(F)
            @test ninf_code(F) == UInt8(2^K - 1)
        else
            @test ninf_code(F) === nothing
        end
    end

    # decode(code_one) == 1.0 for all 120 formats
    for F in FORMAT_ALIASES
        @test decode(fromcode(F, one_code(F))) == 1.0
    end

    # Binary8p1uf range: 2^-127 … 2^126
    @test decode(minpositive(Binary8p1uf)) == 2.0^-127
    @test decode(maxfinite(Binary8p1uf)) == 2.0^126

    # every datum is exactly representable in Float32 and Float64
    for F in FORMAT_ALIASES, c in 0:(2^kbits(F) - 1)
        v = decode(fromcode(F, c))
        isfinite(v) || continue
        @test Float64(Float32(v)) == v
    end

    # draft revision pinned and visible
    @test occursin("P3109", DRAFT_REVISION)
end
