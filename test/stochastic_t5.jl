# T5 — stochastic execution (§17.1): entropy golden vectors, determinism,
# entropy contract errors, small-N grids (the projection-level R grids run
# in T2).

using FloatBytes: entropy_word, IndexedEntropy, SLOT_ELEMENT, SLOT_SCALE,
    codepoint, fromcode

@testset "T5 stochastic" begin
    # golden vectors: exact chaining frozen at Stage 8 (§12.1)
    ie = IndexedEntropy(UInt64(0x123456789abcdef0), UInt64(7), UInt64(42))
    @test entropy_word(ie, 0, SLOT_ELEMENT, Val(8)) == 0x0000000000000029
    @test entropy_word(ie, 1, SLOT_ELEMENT, Val(8)) == 0x00000000000000e9
    @test entropy_word(ie, 0, SLOT_SCALE, Val(8)) == 0x000000000000007b
    @test entropy_word(ie, 12345, SLOT_ELEMENT, Val(32)) == 0x00000000ba03b4ca
    @test entropy_word(ie, 3, SLOT_ELEMENT, Val(1)) == 0x0000000000000001
    @test entropy_word(IndexedEntropy(1), 0, SLOT_ELEMENT, Val(16)) == 0x00000000000006b3

    # R < 2^N always; deterministic replay
    for N in (1, 5, 32, 63), i in (0, 1, 99)
        w = entropy_word(ie, i, SLOT_ELEMENT, Val(N))
        @test w < UInt64(1) << N
        @test w == entropy_word(ie, i, SLOT_ELEMENT, Val(N))
    end
    # distinct keys give distinct words (overwhelmingly)
    words = [entropy_word(ie, i, SLOT_ELEMENT, Val(32)) for i in 0:999]
    @test length(unique(words)) == 1000

    # entropy contract: stochastic ρ requires entropy at the seam
    F = Binary8p4se
    x = F(1.0); y = F(0.0625)
    ρs = Projection(StochasticA(3), SatNone())
    @test_throws EntropyRequiredError apply(Add, F, ρs, x, y)
    # N > 32 carriers are a deferred stage: explicit error, never a guess
    @test_throws ArgumentError apply(Add, F, Projection(StochasticA(40), SatNone()),
                                     x, y, entropy = UInt64(0))
    @test_throws ArgumentError StochasticA(0)
    @test_throws ArgumentError StochasticA(64)

    # IndexedEntropy at the seam: reproducible, index-sensitive
    a = apply(Add, F, ρs, x, y, entropy = ie, index = 5)
    b = apply(Add, F, ρs, x, y, entropy = ie, index = 5)
    @test codepoint(a) == codepoint(b)

    # stochastic mean drifts toward the true value: StochasticA over all R
    # rounds x up with probability equal to the fraction position.
    let N = 4, v = 1.0 + 5 / 16 * 0.125            # fraction 5/16 of an ulp (P=4 ⇒ ulp(1)=1/16... use exact)
        F2 = Binary5p2se                            # P = 2, B = 4; ulp at 1 = 1/2
        vv = 1.0 + 5 / 16 * 0.5                    # ν = 5/16
        ups = 0
        for R in 0:15
            r = apply(Convert, F2, Projection(StochasticA(4), SatNone()),
                      quantize(Binary8p6se, vv), entropy = UInt64(R) << 60)
            ups += (FloatBytes.decode(r) > 1.0)
        end
        @test ups == 5                              # ⌊ν·2^N⌋ = 5 of 16 R values round away
    end
end
