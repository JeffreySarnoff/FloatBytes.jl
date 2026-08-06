# execution/stochastic.jl — indexed entropy (§12.1).
# Kernels never consult global RNG state; the entropy word is data, keyed
# by (version, seed, stream, invocation, logical index, projection slot).
# Same key ⇒ same word, on any schedule, thread count, or storage layout.

"""
    IndexedEntropy(seed, stream, invocation)

Reproducible entropy source for stochastic rounding. `entropy_word`
derives one word per (logical index, projection slot) key.
"""
struct IndexedEntropy
    seed::UInt64
    stream::UInt64
    invocation::UInt64
end
IndexedEntropy(seed::Integer) = IndexedEntropy(UInt64(seed), UInt64(0), UInt64(0))

const ENTROPY_ALGO_VERSION = UInt64(0x01)

# Projection slots (§12.1). The catalogue's authoritative slot list; scalar
# operations use SLOT_ELEMENT.
const SLOT_SCALE     = UInt32(0)
const SLOT_ELEMENT   = UInt32(1)
const SLOT_TIE_BREAK = UInt32(2)

# splitmix64 finalizer; canonical constants.
@inline function mix64(z::UInt64)
    z = xor(z, z >> 30) * 0xbf58476d1ce4e5b9
    z = xor(z, z >> 27) * 0x94d049bb133111eb
    return xor(z, z >> 31)
end

"""
    entropy_word(ie, index, slot, ::Val{N}) -> UInt64

The top-`N`-bit entropy value `R < 2^N` for one logical output. `N` is a
`Val` so `64 - N` folds to a constant shift.
"""
@inline function entropy_word(ie::IndexedEntropy, index::Int, slot::UInt32,
                              ::Val{N})::UInt64 where {N}
    h = mix64(ie.seed + 0x9e3779b97f4a7c15 * ENTROPY_ALGO_VERSION)
    h = mix64(xor(h, ie.stream))
    h = mix64(xor(h, ie.invocation))
    h = mix64(xor(h, (UInt64(index) << 32) | UInt64(slot)))
    return h >> (64 - N)
end
