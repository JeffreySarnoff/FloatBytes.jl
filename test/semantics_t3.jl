# T3 — operation kernels vs the reference (§17.1): exhaustive same-format
# binaries for small K, sampled 8-bit spot checks, special-row pins,
# FAA associativity.

using FloatBytes: kbits, fromcode, codepoint, decode, defined_code, FORMAT_ALIASES
using Random: Xoshiro

const T3_RNG = Xoshiro(0x03109)

function ref_op_code(f, op::Symbol, rmode, rsat, codes...)
    vals = [DR.decode_ref(f, c) for c in codes]
    X = op === :add ? DR.omega_add(vals...) :
        op === :sub ? DR.omega_subtract(vals...) :
        op === :mul ? DR.omega_multiply(vals...) :
        op === :div ? DR.omega_divide(vals...) :
        op === :fma ? DR.omega_fma(vals...) :
        op === :faa ? DR.omega_faa(vals...) :
        op === :sqrt ? DR.omega_sqrt(vals...) :
        op === :recip ? DR.omega_recip(vals...) :
        op === :rsqrt ? DR.omega_rsqrt(vals...) :
        op === :hypot ? DR.omega_hypot(vals...) :
        op === :convert ? DR.omega_convert(vals...) :
        error("unknown $op")
    return DR.project_ref(f, rmode, rsat, X)
end

const BIN_OPS = [(Add, :add), (Subtract, :sub), (Multiply, :mul), (Divide, :div),
                 (Hypot, :hypot)]
const UN_OPS = [(Sqrt, :sqrt), (Recip, :recip), (RSqrt, :rsqrt), (Convert, :convert)]

function check_bin(F, f, op, rop, ρ, rmode, rsat, cx, cy)
    x = fromcode(F, cx); y = fromcode(F, cy)
    got = defined_code(op, F, ρ, UInt64(0), x, y)
    want = ref_op_code(f, rop, rmode, rsat, cx, cy)
    got == UInt8(want) || @info "T3 mismatch" F rop cx cy rmode rsat got want
    @test got == UInt8(want)
end

@testset "T3 semantics" begin
    small = [F for F in FORMAT_ALIASES if kbits(F) <= (EXHAUSTIVE ? 5 : 4)]
    rhos = [(Projection(NearestTiesToEven(), SatNone()), :nte, :satnone),
            (Projection(TowardNegative(), SatFinite()), :tn, :satfinite),
            (Projection(TowardPositive(), SatNone()), :tp, :satnone)]

    @testset "binary exhaustive (small K)" begin
        for F in small
            f = rformat(F)
            n = 2^kbits(F)
            for (ρ, rm, rs) in rhos[1:2]
                for (op, rop) in BIN_OPS
                    for cx in 0:n-1, cy in 0:n-1
                        check_bin(F, f, op, rop, ρ, rm, rs, cx, cy)
                    end
                end
            end
        end
    end

    @testset "unary exhaustive (K ≤ 5)" begin
        for F in [F for F in FORMAT_ALIASES if kbits(F) <= 5]
            f = rformat(F)
            for (ρ, rm, rs) in rhos, (op, rop) in UN_OPS
                for c in 0:2^kbits(F)-1
                    x = fromcode(F, c)
                    got = defined_code(op, F, ρ, UInt64(0), x)
                    want = ref_op_code(f, rop, rm, rs, c)
                    @test got == UInt8(want)
                end
            end
        end
    end

    @testset "8-bit spot (sampled)" begin
        for F in (Binary8p4se, Binary8p1uf, Binary8p3sf, Binary8p5ue)
            f = rformat(F)
            for (ρ, rm, rs) in rhos[1:2]
                for (op, rop) in BIN_OPS
                    for _ in 1:(EXHAUSTIVE ? 4000 : 400)
                        cx, cy = rand(T3_RNG, 0:255), rand(T3_RNG, 0:255)
                        check_bin(F, f, op, rop, ρ, rm, rs, cx, cy)
                    end
                end
                for (op, rop) in UN_OPS, c in 0:255
                    x = fromcode(F, c)
                    @test defined_code(op, F, ρ, UInt64(0), x) ==
                          UInt8(ref_op_code(f, rop, rm, rs, c))
                end
            end
        end
    end

    @testset "FMA/FAA (sampled + specials)" begin
        for F in (Binary4p2se, Binary4p3uf, Binary8p4se, Binary8p2ue)
            f = rformat(F)
            n = 2^kbits(F)
            ρ, rm, rs = rhos[1]
            specials = unique([0, 1, n - 1, n - 2, n ÷ 2, n ÷ 2 - 1, n ÷ 2 + 1])
            triples = [(rand(T3_RNG, 0:n-1), rand(T3_RNG, 0:n-1), rand(T3_RNG, 0:n-1))
                       for _ in 1:(EXHAUSTIVE ? 3000 : 500)]
            append!(triples, [(a, b, c) for a in specials for b in specials for c in specials])
            for (ca, cb, cc) in triples
                x, y, z = fromcode(F, ca), fromcode(F, cb), fromcode(F, cc)
                @test defined_code(FMA, F, ρ, UInt64(0), x, y, z) ==
                      UInt8(ref_op_code(f, :fma, rm, rs, ca, cb, cc))
                gotfaa = defined_code(FAA, F, ρ, UInt64(0), x, y, z)
                @test gotfaa == UInt8(ref_op_code(f, :faa, rm, rs, ca, cb, cc))
                # associativity/commutativity of FAA (§4.4, draft §4.10.7)
                @test gotfaa == defined_code(FAA, F, ρ, UInt64(0), z, y, x)
                @test gotfaa == defined_code(FAA, F, ρ, UInt64(0), y, z, x)
            end
        end
    end

    @testset "cross-format Convert" begin
        for Fx in (Binary8p4se, Binary3p1sf, Binary5p2ue, Binary8p1uf),
            Fr in (Binary8p4se, Binary4p2sf, Binary6p3ue, Binary3p2se)
            fr = rformat(Fr)
            fx = rformat(Fx)
            for (ρ, rm, rs) in rhos, c in 0:2^kbits(Fx)-1
                x = fromcode(Fx, c)
                v = DR.decode_ref(fx, c)
                want = DR.project_ref(fr, rm, rs, DR.omega_convert(v))
                @test defined_code(Convert, Fr, ρ, UInt64(0), x) == UInt8(want)
            end
        end
    end

    @testset "special-row pins (draft-vs-IEEE754 differences)" begin
        F = Binary8p4se
        ρ = Projection(NearestTiesToEven(), SatNone())
        nan = fromcode(F, FloatBytes.nan_code(F))
        pinf = fromcode(F, FloatBytes.pinf_code(F))
        ninf = fromcode(F, FloatBytes.ninf_code(F))
        one_ = one(F); zero_ = zero(F)
        @test isnan(one_ / zero_)                       # x/0 = NaN, not Inf
        @test isnan(inv(zero_))                         # Recip(0) = NaN
        @test iszero(inv(pinf))                         # Recip(±∞) = 0
        @test isnan(apply(RSqrt, F, ρ, zero_))          # RSqrt(0) = NaN
        @test iszero(one_ / pinf)                       # x/±∞ = 0
        @test isnan(pinf + ninf)
        @test isnan(apply(FMA, F, ρ, zero_, pinf, one_))  # 0×∞ NaN
        @test isnan(hypot(nan, pinf))                   # NaN rows precede Inf rows
        @test isnan(sqrt(-one_))
        @test decode(apply(Tanh, F, ρ, pinf)) == 1.0
        @test isnan(nextgreaterthan(pinf))              # unlike IEEE nextUp
        # negate on unsigned: NaN unless zero (SatNone)
        U = Binary8p3uf
        @test isnan(-(one(U)))
        @test iszero(-(zero(U)))
        # Minimum/Maximum propagate NaN; Number variants do not
        @test isnan(min(nan, one_))
        @test decode(apply(MinimumNumber, F, ρ, nan, one_)) == 1.0
        @test decode(apply(MaximumNumber, F, ρ, nan, one_)) == 1.0
    end

    @testset "selection family vs decoded semantics" begin
        for F in (Binary5p2se, Binary4p3uf)
            ρ = Projection(NearestTiesToEven(), SatNone())
            n = 2^kbits(F)
            for cx in 0:n-1, cy in 0:n-1
                x, y = fromcode(F, cx), fromcode(F, cy)
                dx, dy = decode(x), decode(y)
                # Minimum
                want = isnan(dx) || isnan(dy) ? NaN : min(dx, dy)
                @test isequal(decode(apply(Minimum, F, ρ, x, y)), want)
                want = isnan(dx) || isnan(dy) ? NaN : max(dx, dy)
                @test isequal(decode(apply(Maximum, F, ρ, x, y)), want)
                # MinimumMagnitude per draft tie rule
                if !isnan(dx) && !isnan(dy)
                    am = abs(dx) != abs(dy) ? (abs(dx) < abs(dy) ? dx : dy) :
                         (dx < dy ? dx : dy)
                    @test decode(apply(MinimumMagnitude, F, ρ, x, y)) == am ||
                          (isinf(dx) && isinf(dy))
                end
                # comparisons match Float64 semantics exactly
                @test (x < y) == (dx < dy)
                @test (x == y) == (dx == dy)
                @test compare_greater_equal(x, y) == (dx >= dy)
            end
        end
    end
end
