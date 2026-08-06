```@meta
CurrentModule = FloatBytes
DocTestSetup = :(using FloatBytes)
```

# Stochastic rounding

The three draft stochastic modes ([`StochasticA`](@ref),
[`StochasticB`](@ref), [`StochasticC`](@ref)) consume `N` random bits
per rounding decision. **Entropy is input**: kernels never touch global
RNG state, so results are reproducible and exhaustively testable by
construction. Requesting a stochastic ρ without an entropy source throws
[`EntropyRequiredError`](@ref) up front.

Two ways to supply entropy to [`apply`](@ref) / [`quantize`](@ref):

- a raw `UInt64` word — its **top `N` bits** are the entropy value `R`;
- an [`IndexedEntropy`](@ref), which derives the word from the key
  `(seed, stream, invocation, index, slot)` via a splitmix64 chain.

```jldoctest stoch
julia> using FloatBytes

julia> ρ = Projection(StochasticA(8), SatNone());

julia> ie = IndexedEntropy(42);

julia> r = quantize(Binary5p2ue, 1.06, ρ; entropy = ie, index = 0)
Binary5p2ue(0x10 ↦ 1.0)

julia> r == quantize(Binary5p2ue, 1.06, ρ; entropy = ie, index = 0)  # replayable
true
```

Indexed entropy is what makes threaded stochastic execution
schedule-independent: the word depends only on the logical key, never on
evaluation order. An RNG handed into kernels could not have that
property.

## The averaging law

`StochasticA(N)` rounds away with probability `⌊ν·2^N⌋ / 2^N` where `ν`
is the discarded fraction — over all `2^N` entropy values the expected
result is (to within the floor) the exact value:

```jldoctest stoch
julia> v = Binary8p6se(1.15625);                # ν = 5/16 of a Binary5p2se ulp

julia> count(R -> Float64(apply(Convert, Binary5p2se,
                       Projection(StochasticA(4), SatNone()), v,
                       entropy = UInt64(R) << 60)) > 1.0, 0:15)
5
```

## Contract details

- Deterministic modes ignore the entropy word; stochastic modes consume
  exactly one word per logical output, **including for exact and special
  results** — fixed consumption is what makes slicing and replay
  auditable.
- `entropy = 0` under a stochastic mode is a legal value, not a
  sentinel.
- `N ∈ 1:63` is validated at mode construction. This build implements
  `N ≤ 32` end to end; `N > 32` throws an `ArgumentError` at the seam
  (the exact-residue carriers it requires are a deferred implementation
  stage — an explicit error, never a silent approximation).
