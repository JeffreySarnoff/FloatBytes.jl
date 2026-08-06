# Grep-enforced structural rules (§17.6, the subset applicable to the
# implemented stages).

const SRC = joinpath(@__DIR__, "..", "src")
const REF = joinpath(@__DIR__, "..", "reference")

function srcfiles(root)
    out = String[]
    for (dir, _, files) in walkdir(root), f in files
        endswith(f, ".jl") && push!(out, joinpath(dir, f))
    end
    return out
end

@testset "structural rules" begin
    # reference/ never imports src/
    for f in srcfiles(REF)
        s = read(f, String)
        @test !occursin("using FloatBytes", s)
        @test !occursin("import FloatBytes", s)
    end

    # MPFR machinery (directed rounding / precision control) under src/
    # only in operations/elementary.jl
    for f in srcfiles(SRC)
        endswith(f, joinpath("operations", "elementary.jl")) && continue
        s = read(f, String)
        @test !occursin("setrounding", s)
        @test !occursin("setprecision", s)
    end

    # reinterpret into a FloatByte only in domain.jl (and the future dense
    # lowering); external/adapters reinterpret host floats, never FloatByte
    for f in srcfiles(SRC)
        endswith(f, "domain.jl") && continue
        s = read(f, String)
        @test !occursin(r"reinterpret\((T|F|Fr|typeof)", s)
    end

    # the string SmallFloats never appears under src/
    for f in srcfiles(SRC)
        @test !occursin("SmallFloats", read(f, String))
    end

    # load-order comment in FloatBytes.jl matches the include order
    mod = read(joinpath(SRC, "FloatBytes.jl"), String)
    includes = [m.captures[1] for m in eachmatch(r"include\(\"([^\"]+)\"\)", mod)]
    header = split(mod, "module FloatBytes")[1]
    for inc in includes
        @test occursin(basename(inc), header)
    end

    # docstring escape rule (§5.6): no \( \[ \text in src docstrings
    for f in srcfiles(SRC)
        s = read(f, String)
        @test !occursin("\\text", s)
    end
end
