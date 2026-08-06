```@meta
CurrentModule = FloatBytes
```

# API reference

## Index

```@index
Pages = ["reference.md"]
```

## Types and formats

```@docs
FloatByte
SignedFormat
UnsignedFormat
FiniteDomain
ExtendedDomain
FORMAT_ALIASES
```

## Construction and code identity

```@docs
quantize
fromcode
codepoint
isvalidcode
decode
```

## Projection

```@docs
Projection
DEFAULT_PROJECTION
with_projection
NearestTiesToEven
NearestTiesToAway
TowardPositive
TowardNegative
TowardZero
ToOdd
StochasticA
StochasticB
StochasticC
SatFinite
SatPropagate
SatNone
```

## Applying operations

```@docs
apply
defined_code
```

## Operation singletons

Each singleton names its draft clause; all are used through
[`apply`](@ref).

```@autodocs
Modules = [FloatBytes]
Filter = t -> t isa FloatBytes.Operation
```

## Format-level and value queries

```@docs
bitwidth
signedness
domain
exponent_bitwidth
trailing_significand_bits
exponent_bias
maxfinite
minfinite
minpositive
maxsubnormal
minnormal
```

## Comparison, classification, neighbours

```@docs
compare_less
compare_less_equal
compare_equal
compare_greater_equal
compare_greater
total_order
floatclass
FloatByteClass
nextgreaterthan
nextlessthan
```

```@autodocs
Modules = [FloatBytes]
Filter = t -> t isa FloatBytes.FloatByteClass
```

## Stochastic execution

```@docs
IndexedEntropy
entropy_word
```

## Errors

```@docs
InconclusiveProjectionError
EntropyRequiredError
WorkspaceError
AliasingError
```

## Format aliases

The 120 aliases follow one naming rule (`Binary{K}p{P}{s|u}{e|f}`) and
are enumerated by [`FORMAT_ALIASES`](@ref); each carries a docstring
with its parameters (try `?Binary8p4se` in the REPL).

```@autodocs
Modules = [FloatBytes]
Filter = t -> t isa Type && t <: FloatBytes.FloatByte && isconcretetype(t)
```
