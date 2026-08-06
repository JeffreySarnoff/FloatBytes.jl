# show.jl — display discipline (§5.6). `show` never throws and never
# projects: alias, code point, decoded value from the atlas. Out-of-range
# bytes fabricated by raw `reinterpret` are masked and flagged, not errors.

function _alias_name(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D}
    string("Binary", K, "p", P, S === SignedFormat ? "s" : "u",
           D === ExtendedDomain ? "e" : "f")
end

function Base.show(io::IO, x::FloatByte{K,P,S,D}) where {K,P,S,D}
    T = FloatByte{K,P,S,D}
    c = codepoint(x)
    name = _alias_name(T)
    if !isvalidcode(T, c)
        masked = UInt8(c & ((1 << K) - 1))
        print(io, name, "(0x", string(c, base = 16, pad = 2), " INVALID ↦ masked 0x",
              string(masked, base = 16, pad = 2), ")")
        return
    end
    v = decode(x)
    print(io, name, "(0x", string(c, base = 16, pad = 2), " ↦ ", v, ")")
end

Base.show(io::IO, ::MIME"text/plain", x::FloatByte) = show(io, x)
