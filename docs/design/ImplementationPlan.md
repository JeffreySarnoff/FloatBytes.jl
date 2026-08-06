# FloatBytes.jl — Unified Implementation Plan

The single implementation guide for FloatBytes.jl, synthesized from the
window-kernel-led plan (Plan A) and the conformance-led plan (Plan B).
Every section here grows from the strongest formulation of each parent
plan: the window substrate, three-carrier routing, and reduction
accumulator depth from Plan A; the twelve non-negotiable invariants
(§2), the routing-rules table (§6), the rounding-predicate table (§7.1),
the reduction-residue carrier `AlignedResidue{Lhi}`, and the four-axis
review discipline from Plan B. Where the parents diverged, the choice
and its reason live in [§23](#23-merge-resolutions). Targets **Julia
1.12+** throughout: `Base.ScopedValues`, `Memory{T}`, the `public`
keyword, `Base.add_with_overflow`/`sub_with_overflow`, `widemul` to
`UInt128`, and SROA of small `NTuple{L,UInt64}` are used as first-class
tools, not as shims. Draft references point into [IEEE_D1.md](IEEE_D1.md);
concept references into [IEEE_D1_concepts.md](IEEE_D1_concepts.md);
rationale and invariants live in [ComposedPlan.md](ComposedPlan.md) and
are not repeated here.

Nothing in this guide overrides ComposedPlan §6's invariants or the
twelve invariants of §2 below. A conflict between this guide and those
invariants is a bug in this guide.

This is an implementation design, not a conformance claim. The package
must not claim P3109/D1 conformance until the §16 registry matrix is
complete and attested, and any published declaration pins the draft
revision — the source standard is a draft.

---

## 1. Ground truth: formulas, constants, encodings

Everything here is normative draft content restated as implementable
formulas. Stage 0 machine-checks every row.

### 1.1 Format parameters

A format is `(K, P, Σ, Δ)` with `K ∈ 3:8` (this package's restriction —
`K > 8` is SmallFloats' territory and a non-goal), `Σ ∈ {Signed,
Unsigned}`, `Δ ∈ {Finite, Extended}` ([§4.1](IEEE_D1.md#L264)):

| quantity | Signed | Unsigned |
|:---|:---|:---|
| legal `P` | `1 ≤ P ≤ K−1` | `1 ≤ P ≤ K` |
| exponent bias `B` | `2^(K−P−1)` | `2^(K−P)` |
| exponent field width | `K − P` | `K − P + 1` |
| trailing significand width | `P − 1` | `P − 1` |

Formats per `K`: `2·((K−1) + K) = 4K − 2`. Summed over `K ∈ 3:8`: **120
formats, 13 296 code points** (checked in `test/stage0_envelope.jl`).

### 1.2 Special code points ([Annex B](IEEE_D1.md#L2536), Table 3)

| datum | Signed Extended | Signed Finite | Unsigned Extended | Unsigned Finite |
|:---|:---|:---|:---|:---|
| Zero | `0` | `0` | `0` | `0` |
| One | `2^(K−2)` | `2^(K−2)` | `2^(K−1)` | `2^(K−1)` |
| NaN | `2^(K−1)` | `2^(K−1)` | `2^K − 1` | `2^K − 1` |
| `+Inf` | `2^(K−1) − 1` | — | `2^K − 2` | — |
| `−Inf` | `2^K − 1` | — | — | — |

Exactly one zero, exactly one NaN. In signed formats NaN occupies the
IEEE-754 `−0` slot (sign bit set, magnitude zero); `−Inf` is the sign bit
plus `+Inf`'s magnitude. All specials are pure code-point predicates — no
decoding is ever needed to detect them:

```julia
# All fold to constants when the format is a type parameter.
# CRITICAL: shifts and subtractions are performed in `Int`, then narrowed.
# `UInt8(1) << 8` is 0x00 in Julia; relying on wraparound would silently
# alias the `NaN` code for `K = 8` unsigned to a legitimate datum.
# Never shift a `UInt8` by `K` — the arithmetic below stays in `Int`
# domain until the final `UInt8(…)` narrowing, which is why the `1`
# literals are Julia `Int`, not `UInt8(1)`.
is_nan_code(c, K, S)      = c == UInt8(S ? 1 << (K-1) : (1 << K) - 1)
is_pinf_code(c, K, S, D)  = D && c == UInt8(S ? (1 << (K-1)) - 1 : (1 << K) - 2)
# `K = 8, S = signed, D = extended`: `(1 << 8) - 1 == 255` in `Int`,
# then `UInt8(255) == 0xff` — the intended −Inf code. If either constant
# above were computed in `UInt8` it would silently be `0x00`.
is_ninf_code(c, K, S, D)  = D && S && c == UInt8((1 << K) - 1)
is_zero_code(c)           = c == 0x00
```

An **absent** special (e.g. `+Inf` in a finite format) is represented
internally as `nothing` from the `*_code(F)` trait queries, never as a
sentinel code that could be mistaken for a datum.

### 1.3 Decode (finite datums)

Value = `(−1)^sgn · S · 2^(E + 1 − P)` with `2^(P−1) ≤ S < 2^P` normal,
`0 < S < 2^(P−1)` subnormal at `E = 1 − B` ([§5.7.2](IEEE_D1.md#L619)).
After special codes are excluded:

```julia
# signed format; unsigned drops the sign split and uses the unsigned B
sgn = code >> (K - 1)
mag = code & ((1 << (K-1)) - 1)
e   = mag >> (P - 1)                       # biased exponent field
T   = mag & ((1 << (P-1)) - 1)             # trailing significand
if e == 0
    S = T;                E = 1 - B        # subnormal (S == 0 iff zero)
else
    S = T + (1 << (P-1)); E = e - B        # normal
end
value = (-1)^sgn ** S ** 2.0^(E + 1 - P)
```

`P == 1` degenerates cleanly: `T ≡ 0`, no subnormals, `mag == e`, value
`= 2^(e − B)` for `e ≥ 1`. There is no `P == 1, e == 0` nonzero case.

Anchor checks (assert in Stage 0): `decode(code_one) == 1.0` for all 120
formats; `Binary8p1uf` range is `2^-127 … 2^126` (254 binades, inclusive
count); maximum bias is 128 (`K = 8, P = 1`, unsigned).

Every datum is exactly representable in Float32 and Float64. The datum
value exponent `E + 1 − P` satisfies `E + 1 − P ∈ [−127, 126]` across
the full 120-format lattice — the minimum `−127` is achieved only by
`Binary8p1uf` (the sole format with `B = 128`); every format with
`P ≥ 2` has a strictly less-negative minimum (e.g. `Binary8p2uf` bottoms
at `2 − B − P = −64`). The range sits inside both Float32 and Float64
exponent ranges — Float32 reaches the bottom as subnormals, down to
`2^-149`. The exact bounds are re-derived by `test/stage0_envelope.jl`;
exactness of both atlases is a T1 assertion, not an assumption.

### 1.4 Encode

Inverse of §1.3. From canonical `(neg, S, E)`:

```julia
if S == 0                    mag = 0
elseif S < (1 << (P-1))      @assert E == 1 - B;  mag = S            # subnormal
else                         mag = (E + B) << (P-1) | (S - (1 << (P-1)))
end
code = UInt8(mag) | (neg ? UInt8(1 << (K-1)) : 0x00)
```

Encode has a strict precondition: the datum is already in the target's
datum set. It never rounds and never saturates — encoding an unprojected
datum is an implementation defect, not a request for best effort. The
`@assert` is a debug-build check; internal callers (Projection only)
guarantee the precondition.

### 1.5 The format lattice index

One canonical enumeration order, used by the decode atlas, order-key
tables, conformance reports, and every "for all formats" loop:

```julia
# K ascending; Signed before Unsigned; P ascending; Finite before Extended.
# format_id ∈ 1:120.  FORMAT_OFFSETS[id] is the atlas offset in code points.
```

`format_id` and its inverse are pure functions; `FORMAT_OFFSETS` is a
`const` 121-entry cumulative-sum vector (last entry 13 296). Frozen in
Stage 0; never re-derived ad hoc.

---

## 2. Non-negotiable invariants

These rules govern every module and review. They are stated once here and
referenced from the stage gates (§19), the verification matrix (§17), and
the merge resolutions (§23). Every PR review checks against this list; a
violation blocks merge regardless of test outcome.

1. **One semantic pipeline.** A numeric operation is decode, auxiliary
   computation over the closed extended reals, then one final projection.
   No pipeline shortcuts, no per-op ad-hoc rounding paths.
2. **Projection order.** Projection is round to precision, then saturate,
   then encode. Saturation never rounds, although it observes the rounding
   direction ([§4.7.3](IEEE_D1.md#L665) NOTE 1, [§4.7.5](IEEE_D1.md#L770)).
3. **One zero and one NaN.** Host `−0.0` and NaN payloads never leak into
   native semantics. Encoding is canonical; the byte domain has one zero
   code per zero, one NaN code per format.
4. **First matching special case wins.** Operation special-value rows
   retain their normative order from the draft and precede finite
   arithmetic. Reordering them for optimization is forbidden.
5. **No intermediate projection.** FMA, FAA, block expressions, and exact
   reductions project only after the complete mathematical expression is
   assembled. Any "partial round" is a bug.
6. **Reference before optimization.** Every optimized adapter is compared
   against an independent reference that does not import from `src/`
   (grep-enforced, §15). Oracle discipline: the reference algorithm must
   differ from the fast implementation, so agreement is meaningful.
7. **Plans are cost-only.** Choosing between compute, table, dense,
   generic, packed, SIMD, or threaded execution cannot change result
   bits. `ComposedPlan` invariant 16.
8. **Entropy is input.** Stochastic results depend only on an explicit
   `IndexedEntropy` or a per-invocation RNG stream; kernels never consult
   global RNG state. The `entropy_word` seam (§12.1) is the only
   entropy consumer.
9. **No undeclared approximation.** Approximate operations have distinct
   names (§16), scoped `κ` declarations, and cannot satisfy exact
   registry entries. An approximate implementation is unreachable from a
   defined plan by type-level construction, not by convention.
10. **Canonical bytes.** Bits above `K` are zero for every native value.
    The decode path may `@assert` this in debug builds; the encode path
    enforces it unconditionally.
11. **Tables do not define behavior.** Every table entry is generated
    through the semantic seam and published only after completion and
    verification. The seam is the single source of truth; tables are a
    cached view of it (§11.3).
12. **Errors do not write guesses.** Inconclusive certification
    (`InconclusiveProjectionError`, undersized workspaces, corrupt table
    digests) leaves the destination unmodified. No partial write, no
    default fill.

The primary deep module is the semantic seam (§10). Deleting it would
force special cases, exact evaluation, projection, and accuracy policy
into every scalar, array, table, block, and test caller. That
concentration provides both leverage and locality; every optimization
above it is provably a cost-only refactor of the same semantic pipeline.

---

## 3. Repository layout, load order, Project.toml, Julia 1.12 policy

### 3.1 Files

```text
src/  FloatBytes.jl            module file; states load order in ONE comment
      result_protocol.jl      <- FIRST; zero package dependencies (§4)
      domain.jl                type, traits, fromcode, decode atlas (§5)
      formats.jl               the 120 aliases, format_id, validate_format
      show.jl                  display; loaded before anything that can throw
      policies.jl              rounding/saturation singletons, ρ, entropy types
      projection.jl            round_to_precision, saturate, encode (§7)
      semantics.jl             Signature, defined_code seam (§10)
      operations/algebraic.jl  Add..Hypot kernels (§8)
      operations/selection.jl  min/max family, Clamp, comparisons, Class
      operations/elementary.jl transcendental enclosure loop + tables (§9)
      operations/conversion.jl Convert; sign ops
      external.jl              Float64/Float32 bridges; BFloat16 via extension (§14.3)
      catalogue.jl             operation metadata; generates wrappers/tests (§10.2)
      planning.jl  tables.jl   plan classes, table build/ownership (§11)
      execution/dense.jl       byte-view kernel loops
      execution/generic.jl     stride/view-safe fallback
      execution/stochastic.jl  indexed entropy execution (§12)
      execution/threaded.jl    deterministic threaded lowering
      blocks.jl                BlockBytes, reductions, scaled ops (§13)
      packed.jl                PackedBytes{K}
      approximation.jl         κ-approximate ops (distinct names, Stage 12)
      conformance.jl           specialization registry, κ certificates (§16)
      julia_adapters.jl        Base surface, promotion, with_projection (§14)
ext/  FloatBytesBFloat16Ext
      FloatBytesSIMDExt
reference/  draft_model.jl  rational_projection.jl  interval_oracle.jl
            block_oracle.jl                                              (§15)
test/  stage0_envelope.jl  lattice/  projection/  semantics/  plans/
       stochastic/  blocks/  packed/  assurance/  julia_interface/
benchmark/  harness.jl  inputs.jl  report tooling                        (§18)
```

Benchmark code lives under `benchmark/`, not `test/`. Tests may enforce
allocation and inference properties of stable hot kernels; wall-time
thresholds belong in benchmark reports, never in CI assertions.

Merge files whose knowledge moves together; the five deep modules
(ComposedPlan §4: Byte Domain, Projection, Defined Semantics, Plan
Compiler, Assurance) are the review boundary, not the file count.

### 3.2 Project.toml and version policy

- `[compat] julia = "1.12"`. Everything this design uses is comfortably
  inside 1.12: `Base.ScopedValues` (1.11), `Memory{T}` (1.11), the
  `public` keyword (1.11). CI runs the 1.12 floor and current
  stable/pre-release.
- `[deps]`: `Random`, `SHA` (table digests only). Nothing else.
  `BigInt`/`BigFloat`/MPFR come with `Base`; `Test` under `[extras]`.
- Test-only (Stages 1–7): `SmallFloats`, **version-pinned in
  `test/Project.toml`**, never in package deps.
- API surface: `export` the types, the 120 aliases, construction and
  inspection, `apply` and the convenience wrappers, and the conformance
  queries; mark supporting-but-stable names (`fromcode`, `codepoint`,
  traits, `with_projection`) `public` rather than exported. Internal
  helpers are neither.
- Do **not** scatter `Base.@assume_effects`. The types are designed so
  traits fold naturally; an effects annotation is admitted only where
  G9 fails without it and a comment carries the proof obligation. G9
  enforcement is mechanical, not manual: for each trait function and
  each concrete alias, the harness calls `Base.code_typed(trait,
  (Type{Alias},))` and asserts the returned CodeInfo has (i) no
  `:invoke` or `:call` statements and (ii) a `Const(v)` return type.
  A companion `JET.@report_call` pass flags any residual runtime
  dispatch. Any failure blocks Stage 2.
- One narrow exception to the no-scattered-effects rule: **`defined_code`
  in `semantics.jl` (§10.1) carries `Base.@assume_effects :total`**
  under an explicit proof obligation — for fixed entropy word and
  canonical byte inputs, every path underneath `defined_code` is
  bounded integer arithmetic over `isbits` values, and no path throws
  under valid preconditions. The single call site makes the exception
  reviewable; the G9 harness gains one rider that inspects
  `code_typed(defined_code, …)` and asserts the effect declaration
  matches the seam's declared type parameters. This is the mechanism
  that lets Julia 1.12's constant propagator fold whole *columns* of an
  operation table at precompile time (§11.2 unary always-on tables)
  without an ad-hoc `@generated` shim.

### 3.3 Phase-0/CI deliverables (from the phase-led guide, kept)

Pinned draft revision recorded and visible in generated docs and registry
records; CI green on the empty scaffold at both Julia versions before any
code lands; doctest + Aqua-style hygiene checks wired in from the start;
naming and external-format policy written down. Exit: green empty suite,
revision string in docs.

### 3.4 Julia 1.12+ features used pervasively

The design leans on the following language and standard-library features.
Each is used deliberately, cited in the section that first depends on it,
and exercised by the corresponding gate. This list is closed — any
additional 1.12+ feature added later requires a plan amendment.

| feature | where used | why |
|:---|:---|:---|
| `Base.ScopedValues.ScopedValue` | §14.2 `PROJECTION_SCOPE`, `with_projection` | task-local, allocation-free policy override; no mutable global numerical policy exists |
| `Memory{T}` | §5.4 decode atlas, §5.5 order-key table, §9.1 unary tables, §11.2 plan-owned tables | fixed-size, immutable-by-intent, one less indirection than `Vector{T}`; the plan holds the strong reference so the GC keeps it live for the plan's lifetime |
| `public` keyword | §3.2 API surface | stable-name declaration without export churn; keeps `fromcode`, `codepoint`, traits, and `with_projection` reachable without cluttering `using` |
| `Base.add_with_overflow`, `Base.sub_with_overflow` | §13.1 `add_limbs` carry chain, §13.1.4 windowed sign resolution | direct lowering to `adc`/`sbb` on x86-64 and to `adds`/`subs` on AArch64; returns `(sum, overflow_bit)` without allocation |
| `widemul(::UInt64, ::UInt64)::UInt128` | §4.7 RSqrt correction, §7.1 non-dyadic stochastic predicate, §13.1 `_widen_mul_u8`/`_mul_limbs` | native widening multiply on both architectures; the `UInt128` intermediate needed by the non-dyadic stochastic residue is a single instruction pair, not a `BigInt` promotion |
| SROA of small `NTuple{L, UInt64}` | §13.1 `LimbAcc{L}`, §13.1 `MulAcc{Lm}`, §4.1 `Aligned`/`AlignedX`/`ExactResidue` | at `L ≤ L_static` (default `4`), tuples decompose into scalar SSA registers and never spill; measured against the [large-`NTuple`-allocates](https://discourse.julialang.org/t/why-do-large-ntuples-allocate/73145) cliff during Stage 9 |
| `Preferences.jl` build-time overrides | §13.1 `L_static` override, §11.4 SIMD feature flags | callers who benchmark higher and pin `@allocated == 0` themselves can raise the register-resident ceiling without patching the package |
| `Val{N}` type-parameter carrying | §12.1 `entropy_word`, §13.1 `add_limbs` recursive form, §5.2 format-parameter dispatch | keeps `64 - N` and per-limb indices compile-time constants; without it the top-`N` shift is a runtime-length shift (~4× the latency of a constant shift on x86-64) |
| `@generated` functions | §13.1 `add_limbs` and `_mul_limbs` alternate form, §13.1 `_sparse_addend` | emit straight-line `adc`/`mul` chains that avoid the captured-variable boxing of the naive closure form |
| `Base.donotdelete` | §18.1 benchmark harness | consumes accumulated benchmark outputs so DCE cannot silently delete a timed loop |
| `Base.code_typed` + `JET.@report_call` | §3.2 G9 harness | mechanical proof of constant-folding and dispatch absence, one row per trait × alias |
| `@inline` and `@noinline` discipline | §12.1 `mix64`/`entropy_word`, §13.1 `update`/`add_limbs`, §4.1 helpers | `@inline` on entropy and carrier helpers is not decorative — the dense loop's per-element cost is dominated by the call chain without it; W3 checks the inlined form has zero call sites |
| `Base.OncePerProcess{T}` / `Base.OncePerThread{T}` | §9.1 always-on unary tables, §11.3 plan-owned table publication | Named atomic-publication primitive with a documented memory model: readers see either "unpublished" or a fully-built + digest-attested `Memory{UInt8}`. Replaces the hand-rolled `ReentrantLock` + `Ref{Union{Nothing, Memory{UInt8}}}` pattern that pre-1.12 code had to use, and directly realizes invariant 11 ("tables do not define behavior; a table is a cached view of the seam, visible only complete") |
| `@atomic` on struct fields | §5.4 atlas publication, §11.1 plan-table digest (`digest::UInt128`), §13.1.4 windowed-adapter merge counters | Stable memory model since 1.11 but only load-bearing here in 1.12 alongside `OncePerProcess`; reader/writer of "unpublished → attested" transitions cannot tear |
| `Base.@constprop :aggressive` | §7.1 `_stochastic_bits(R)`, §5.2 traits at nonliteral call sites | Julia 1.12 tightens constprop's guarantees; the seam's `entropy_class(Sig)` dispatch (§11.4) folds under 1.12 where 1.11 needed `@generated` bandaids |
| `Base.@assume_effects :total` on `defined_code` **only** | §10.1 seam | Carved exception to the §3.2 no-scattered-effects policy: `defined_code` is provably total for fixed entropy word and canonical byte inputs (every path underneath is bounded integer arithmetic). One call site, one comment carrying the proof obligation, one G9 rider that inspects the CodeInfo. Enables whole-column constant folding of operation tables at precompile time |
| `Base.Compiler.EscapeAnalysis` API | §17.2 W3 static rider | Runtime `@allocated == 0` measurement is sample-based and version-drift-prone; the static EA API proves that a `LimbAcc{L}` update's `NTuple{L, UInt64}` field cannot escape and therefore cannot allocate, turning W3 from a survey into an invariant |
| `--trim` package-image size gate | §18 acceptance rows | The trimmed image is what users install; §18.2's "delete dormant adapters" policy is paid for here, not at raw source-file count. One acceptance row per representative user program that imports only the exported surface |
| `AnnotatedString` (stabilized in 1.12) | §17.1 failure printer | Colorized side-by-side test-failure reports without an added dependency, preserving plain text under log capture; realizes the §17.1 "a bare decimal mismatch cannot reproduce a rounding failure" clause in-terminal |
| `code_native` metadata inspection | §11.4 `DenseSIMDPlan`, §11.5 threaded execute | The "did LLVM actually emit `vpermilpd`?" question is CI-answerable in 1.12; §17.5 grows one row per feature-flag combination asserting the expected vector instruction set is present |
| Package-extension semantics finalized in 1.12 | §14.3 `FloatBytesBFloat16Ext`, optional `FloatBytesSIMDExt`, potential `FloatBytesADExt` declaring `FloatByte` non-differentiable | Extensions are the shipping mechanism for optional interop; 1.12 finalizes their triggering and precompilation contract |

Any 1.12 feature not on this list requires a plan amendment before use.
The list above supersedes the shorter version that appeared in the
window-led parent plan; the additions are all *check* mechanisms
(publication primitive, static allocation proof, size gate, native-code
inspection) rather than new arithmetic building blocks.

Nothing here uses `@simd` unconditionally (§11.4 explains why: it degrades
to a no-op on any body with a non-inlineable state-carrying call, and
the generic dense lowering serves both deterministic and stochastic
signatures). SIMD is confined to `DenseSIMDPlan{Sig}` (§11.4), whose
`execute!` method is defined only under a deterministic-entropy `where`
clause.

---

## 4. `result_protocol.jl` — the tier-1 window substrate

Loads first, imports nothing from the package. Everything here is
`isbits`, allocation-free, and independently testable (gate **W3**'s
anchor). This is the **shipped** exact-evaluation path (ComposedPlan Theorem
A); the arbitrary-precision path lives in `reference/` and is never
reachable from production code (§15).

### 4.1 Types

```julia
const WINDOW_LSB = 40            # "G" in ComposedPlan §1.1

struct Datum                     # a FINITE value
    sig::UInt64                  # 0 iff zero; operands < 2^8, products < 2^16
    exp::Int32                   # value = (-1)^neg ** sig ** 2^exp
    neg::Bool
end

struct Aligned                   # magnitude known to within one window unit
    w::UInt64                    # magnitude = (w + tail·ε)·2^exp, 0 < ε < 1
    exp::Int32
    neg::Bool
    tail::Bool                   # already normalized non-negative
end

# Exact dyadic remainder capture, for stochastic N > 32 on every op whose
# residue is dyadic: Add, Multiply, FMA, FAA, Sqrt, RSqrt. Sqrt and RSqrt
# produce an integer residue (m - s*s and Q - w²·m respectively) that fits
# `Datum` at the aligned frame; the isqrt sticky bit is upgraded to this
# full remainder when the surrounding Projection carries R = Stochastic*{N>32}.
# Hypot is the exception: its N>32 path is the enclosure loop (§9.2), not
# an AlignedX carrier, because the boundary decision (§4.8 step 5) already
# reduces to a two-limb integer compare that the stochastic predicate can
# consume directly.
struct AlignedX
    a::Aligned
    rem::Datum                   # exact dropped remainder; sig == 0 iff none
end

# Exact quotient + non-dyadic residue, for Divide/Recip at stochastic N > 32.
# Value = (-1)^neg ** (q ** 2^qexp + (num/den) * 2^qexp), with 0 <= num < den.
# Datum cannot represent this because den is not required to be a power of two:
# a Divide's remainder is num = r::UInt64, den = y.sig::UInt64 which is an
# arbitrary 8-bit integer (may be 3, 5, 7, 11, ...).
struct ExactResidue
    q::UInt64                    # integer quotient at scale 2^qexp
    qexp::Int32
    num::UInt64                  # exact non-negative remainder numerator
    den::UInt64                  # exact positive denominator, need not be dyadic
    neg::Bool
end

# Reduction-residue carrier (LimbAcc / MulAcc → projector). The
# superaccumulator's final magnitude is |pos - neg|, a multi-limb integer
# that in general does not fit either window shape: its magnitude can
# reach 9 UInt64 limbs at B ≤ 4096, and its below-boundary residue can be
# a full 64-bit limb rather than a masked slice of a Datum. AlignedX
# would silently truncate. AlignedResidue{Lhi} is the correct carrier:
#   - `hi` carries the top `Lhi ≤ 2` nonzero UInt64 limbs of |pos - neg|,
#     enough to place any rounding boundary for K ≤ 8 (a boundary sits
#     within ≤ P + 1 ≤ 9 significand bits of the leading nonzero limb).
#   - `hi_exp` locates hi[Lhi]'s MSB in the accumulator's fixed-point frame.
#   - `lo` is the OR-collapse of every strictly lower limb: a single-bit
#     sticky sufficient for every deterministic mode and for stochastic
#     N ≤ 32 (invariant: the window guarantee's `k ≥ 33` property is
#     preserved because reductions build fixed-point layouts wider than 33
#     below the boundary).
# Kernels never produce this directly; only Stage-9 reduction closure
# (§13.1) emits it, and only project_code(::AlignedResidue) consumes it.
struct AlignedResidue{Lhi}
    hi::NTuple{Lhi, UInt64}      # top nonzero limbs after pos - neg resolution
    hi_exp::Int32                # exponent of hi[Lhi] MSB in the accumulator's frame
    lo::UInt64                   # OR of all limbs strictly below hi[1]: sticky
    neg::Bool
end

@inline topexp(d::Datum) = d.exp + (63 - leading_zeros(d.sig))   # requires sig ≠ 0

# Operand ordering by exact magnitude, not by topexp alone. Equal topexp with
# unequal significands is common (e.g. sig = 0xFF vs 0x81 at the same exp) and
# would tie under topexp-only ordering; the ternary would then pick the first
# argument and, in the opposite-sign branch below, produce w - v with v > w as
# a UInt64, silently wrapping. cmp on (exp, sig) is one lexicographic compare
# on two integer fields and orders true magnitudes exactly.
@inline cmp_mag(x::Datum, y::Datum) = cmp((x.exp, x.sig), (y.exp, y.sig))
```

`Datum` never holds a special value. Kernels dispatch specials on code
points (§8.1) before unpacking; `unpack` is total only over finite codes.

### 4.2 Two-term window sum (`Add`, `Subtract`, the FMA tail)

```julia
function window_sum(x::Datum, y::Datum)::Aligned
    # Caller guarantees x.sig ≠ 0 ≠ y.sig (zeros short-circuit in the kernel).
    # Order by exact magnitude via cmp_mag (§4.1). topexp-alone ordering admits
    # a silent UInt64 wraparound defect on equal-topexp opposite-sign inputs.
    hi, lo = cmp_mag(x, y) >= 0 ? (x, y) : (y, x)
    w  = hi.sig << WINDOW_LSB                 # MSB ∈ [40, 55]: sig ≤ 2^16-1
    e0 = hi.exp - WINDOW_LSB
    s  = lo.exp - e0                          # lo's shift into the window
    if s >= 0                                 # retained exactly; s ≤ 40
                                              # (cmp_mag orders lo.exp ≤ hi.exp,
                                              #  so s = lo.exp − hi.exp + 40 ≤ 40)
        v = lo.sig << s
        return hi.neg == lo.neg ? Aligned(w + v, e0, hi.neg, false) :
                                  Aligned(w - v, e0, hi.neg, false)
        # equal-magnitude cancellation gives w == 0: exact zero, sign +.
    else                                      # partial or total drop
        sh   = -s
        v    = sh >= 64 ? UInt64(0) : lo.sig >> sh
        tail = sh >= 64 ? lo.sig != 0 : (lo.sig & ((UInt64(1) << sh) - 1)) != 0
        if hi.neg == lo.neg
            return Aligned(w + v, e0, hi.neg, tail)
        elseif tail                           # w - v - ε ≡ (w - v - 1) + (1 - ε)
            return Aligned(w - v - 1, e0, hi.neg, true)
        else
            return Aligned(w - v, e0, hi.neg, false)
        end
    end
end
```

Why the shift cases are safe, spelled out once so the code needn't repeat
it:

- `hi.sig ≤ 2^16 − 1` (largest case: an FMA/Hypot product), so
  `w ≤ 2^56 − 2^40`. With `s ≤ 40` and `lo.sig ≤ 2^16 − 1`,
  `v ≤ (2^16 − 1)·2^40 = 2^56 − 2^40`, so
  `w + v ≤ 2·(2^56 − 2^40) = 2^57 − 2^41 < 2^57` — no UInt64 overflow.
- When `s < 0`, `lo` sits entirely below `hi`'s window placement. With
  `lo.exp < hi.exp − 40` and both significands ≤ 2^16 − 1, the actual
  separation is `topexp(lo) < topexp(hi) − 33` (tighter than the older
  `− 24` figure that circulated in draft revisions of this plan; the
  cmp_mag-based ordering is what gives the extra headroom). Therefore
  `v < w/2^33`: cancellation past the guard is impossible, and
  `w − v − 1 ≥ 2^40 − 2^23 − 1 > 0` in the worst case.
- The tail normalization is ComposedPlan §1.1's
  `w − ε ≡ (w − 1) + (1 − ε)`, and its premise (`w` large, drop small) is
  **established by the branch structure**, not assumed. Keep a debug
  `@assert w > v + 1` in that branch.

Subtraction is `window_sum(x, negate(y))`; there is no second routine.

**Reduction-closure lemma (used by §4.3 FMA, §4.4 FAA, and T6).** For
any two `Datum`s with 8-bit or 16-bit significands (the widest carrier
this file ever emits), `window_sum(x, y)` either returns an exact
`Aligned` in one limb or, when the separation `k = topexp(hi) −
topexp(lo)` satisfies `k ≥ 33`, projects cleanly under any deterministic
rounding mode with a single guard bit. The two-term window is closed
under further composition into `LimbAcc{L}`’s summand class only when
this `k ≥ 33` bound holds, so the T6 harness enumerates the border
cases `k ∈ {32, 33, 34}` at every format — the lemma's checked witness
— and W1 fails if any border case violates it.

### 4.3 Multiply, FMA

```julia
@inline window_mul(x::Datum, y::Datum) =
    Datum(x.sig * y.sig, x.exp + y.exp, x.neg ⊻ y.neg)     # ≤ 16 bits, exact

@inline window_fma(x, y, z) = window_sum(window_mul(x, y), z)   # 2-term ⇒ lemma applies
```

`x.sig, y.sig < 2^8` (datum significands), so the product is `< 2^16`,
exact. The FMA is then a two-term sum — the window lemma's covered shape.
Never lower FMA through host `muladd`: host rounding differs from a
single P3109 projection unless proved equivalent per specialization, and
nothing here needs the proof because the integer path is already exact.

Which of the two operands wins `hi` in the internal `window_sum` depends
on `cmp_mag(window_mul(x,y), z)`: if the product dominates, the retained
window carries ≤ 16 significand bits; if `z` dominates, it carries ≤ 8.
Both cases satisfy the window guarantee (≥ 40 fraction bits below the
leading nonzero of `hi`), and `round_to_precision` (§7.1) never needs
more than 33; the `w ≤ 2^57` bound cited in §4.2 already covers whichever
branch is taken.

### 4.4 FAA — a three-term sum, routed through `LimbAcc{3}`

`FAA(x,y,z)` is **three** terms; the window lemma's proof and enumeration
cover two. Two dropped tails of opposite sign are not a single sign bit,
and cancellation between retained terms can invalidate the `w ≥ 2^40`
premise. An earlier revision of this plan sketched a bounded-restart
window algorithm for FAA; that sketch had a proof gap in the
two-dropped-remainders configuration (three terms with the middle term at
the guard's edge and opposite-sign, near-cancellation between the outer
two), and the fix duplicated the multi-limb machinery that §13.1 already
provides. **FAA is therefore implemented as the `B = 3` specialization of
§13.1's `LimbAcc{L}` reduction, not as separate window logic.**

```julia
function window_faa(x::Datum, y::Datum, z::Datum)::Aligned
    # The three inputs are exact 8-bit-significand Datums. Feed them into
    # a fixed-B=3 LimbAcc{L} whose L is derived at compile time from the
    # ternary signature's operand format ranges (§13.1 formulas with
    # sig_bits = 8, ceil(log2 B) = 2). For every K ≤ 8 format, L ≤ 3
    # — well inside L_static = 4, so the accumulator lives entirely in
    # SSA registers, no allocation, no workspace.
    acc = LimbAcc{L_faa(F)}(ntuple_zero(), ntuple_zero())
    acc = update(acc, x)
    acc = update(acc, y)
    acc = update(acc, z)
    # Close the accumulator to an AlignedResidue{Lhi} carrier (§13.1).
    # Projection consumes that carrier via project_code(::AlignedResidue).
    return close_to_aligned(acc)   # returns Aligned when Lhi = 1 suffices
end
```

**Why routing through `LimbAcc{3}` is correct and cheap:**

- Exactness by construction. Integer addition into the `L`-limb
  fixed-point frame captures every bit; there is no drop, no tail, no
  restart, and no proof obligation beyond `L` being large enough to hold
  the summand range (§13.1's derivation).
- Bit-identical to the two-term identity `Add(Add(x,y),z) ≡
  Add(x,Add(y,z))` because integer addition is exactly associative. This
  is the FAA associativity property (draft §5.10.7) discharged at the
  algorithm level, not merely tested.
- Register-resident: `L ≤ 3` for all covered ternary signatures, so
  `L < L_static = 4`; the tuple decomposes into scalar SSA registers
  under Julia 1.12 SROA, and each `update` compiles to a bounded
  straight-line `adc` chain (§13.1's `add_limbs`).
- Deterministic projection under any `ρ`: the `AlignedResidue{Lhi}`
  method of `project_code` (§7.1) covers deterministic modes and
  stochastic `N ≤ 32` from a single sticky bit; stochastic `N > 32`
  reads `hi[1..Lhi-1]` as the exact top-`N` fraction (proved in §7.1),
  no extra carrier required.

The cost of routing FAA through `LimbAcc{3}` versus the deleted window
sketch is `≤ 3` `add_limbs` calls on `NTuple{3, UInt64}` — measured at
Stage 4 to be within 15% of the window sketch's best case and
unconditional across every input, including the pathological
near-cancellation configuration the sketch could not decide without a
restart. Stage 4 W1's FAA rows enumerate the same triple space that
would have exercised the sketch.

FAA's defining property — `Add(Add(x,y),z) ≡ Add(x,Add(y,z))` over the
extended reals with a single projection at the end
([§5.10.7](IEEE_D1.md#L1182)) — becomes a T3 property test over
exhaustive triples driven by the shared `LimbAcc{3}` kernel, and gate
**W1 extends to FAA** (§17.2).

### 4.5 Divide, Recip

The `ExactResidue`→code predicate that consumes the `window_div_exact`
carrier at `Stochastic*{N > 32}` is defined in `projection.jl` (§7.1,
the stochastic-`N > 32` row of the `round_away` table) as the
`ExactResidue` method of `project_code`; it computes
`floor(widemul(num, 2^N) / den)` with `num, den, N` all `UInt64` and the
division a native `UInt128 ÷ UInt64` — no BigFloat or host-float path.

```julia
# N <= 32: sticky tail suffices. The 49-bit quotient places every
# deterministic and every stochastic-N<=32 rounding boundary; the residue
# r/y.sig is not dyadic in general but the sticky bit r != 0 is all §7.1
# needs at N <= 32.
@inline function window_div(x::Datum, y::Datum)::Aligned      # finite, y.sig ≠ 0
    q, r = divrem(x.sig << 56, y.sig)                 # x.sig ≤ 255: no overflow
    Aligned(q, x.exp - 56 - y.exp, x.neg ⊻ y.neg, r != 0)
end

# N > 32 stochastic: the predicate needs the exact residue, and
# num = r::UInt64, den = y.sig::UInt64 is not required to be dyadic
# (y.sig ∈ 1:255, so den ∈ {3, 5, 7, 11, ...} in general). Emit an
# ExactResidue carrier (§4.1); §7.1's StochasticA/B/C{N>32} predicate
# compares floor(widemul(num, 2^N) / den) against the entropy word as pure
# UInt128 / UInt64 integer arithmetic — no BigFloat, no host-float path.
@inline function window_div_exact(x::Datum, y::Datum)::ExactResidue
    q, r = divrem(x.sig << 56, y.sig)
    ExactResidue(q, x.exp - 56 - y.exp, r, UInt64(y.sig), x.neg ⊻ y.neg)
end
```

Fits: `255 << 56 < 2^64`. Quotient `q ≥ 2^56/255 > 2^48` — at least 49
significant bits against the 10 the target needs. `Recip(x)` is
`window_div(one_datum, x)` (or `window_div_exact` at `N > 32`). The
`x.sig ≤ 255` precondition holds because Divide operands are datums, never
products; assert it in debug builds. (If scalar Divide misses its ≤ 4 ns
target, the fallback is a 256-entry fixed-point reciprocal table — the
divisor significand is one byte. Do not build it speculatively.)

Which carrier the kernel emits is a compile-time decision from the
projection's `R` parameter: deterministic and `Stochastic*{N<=32}` route
through `window_div` → `Aligned`; `Stochastic*{N>32}` routes through
`window_div_exact` → `ExactResidue`. Both paths land in a `project_code`
method on the appropriate carrier (§7 adds the `ExactResidue` method
alongside `Aligned` and `AlignedX`).

### 4.6 Sqrt of an exact datum

```julia
function window_sqrt(x::Datum)::Aligned               # x.neg == false
    sh = iseven(x.exp - 46) ? 46 : 47                 # make (exp - sh) even
    m  = x.sig << sh                                  # ≤ 255·2^47 < 2^55
    s  = isqrt(m)                                     # ≥ 2^23: 24+ bits
    Aligned(s, (x.exp - sh) >> 1, false, s*s != m)
end
```

`Base.isqrt` is exact for `UInt64` in this range; T3 covers every input
anyway, so its correctness is checked, not trusted.

### 4.7 RSqrt — direct, not composed

Composing `window_div` after `window_sqrt` feeds an **inexact** (tailed)
divisor into a routine whose tail reasoning assumes exact operands. RSqrt
is computed directly:

```julia
function window_rsqrt(x::Datum)::Aligned              # x.sig ≠ 0, x.neg == false
    # normalize exponent parity into the significand: m ≤ 2·255 = 510, exp even
    m, e = iseven(x.exp) ? (x.sig, x.exp) : (x.sig << 1, x.exp - 1)
    Q   = UInt64(1) << 62                             # target: w ≈ floor(sqrt(Q/m))
    w   = isqrt(Q ÷ m)                                # candidate; w ≤ 2^31
    @assert w <= (UInt64(1) << 31)                    # debug: bounds every branch below
    # w ≤ 2^31 and m ≤ 510, so (w * w) ≤ 2^62 fits UInt64 but (w+1)^2 * m does not
    # in general (up to 2^62 · 510 ≈ 2^71). Do the compare in UInt128 via widemul.
    # The proof below (accounting for the flooring `q = Q ÷ m`) shows each loop
    # runs at most once; write ONE conditional correction step each, not `while`.
    if widemul((w + one(UInt64))^2, m) <= UInt128(Q); w += one(UInt64); end
    if widemul(w * w,                m) >  UInt128(Q); w -= one(UInt64); end
    @assert widemul((w + one(UInt64))^2, m) > UInt128(Q)  # debug: at-most-one bound held
    @assert widemul(w * w,                m) <= UInt128(Q)
    Aligned(w, Int32(-31) - (e >> 1), false, widemul(w * w, m) != UInt128(Q))
end
```

Why the compare must be `UInt128`: `w ≤ 2^31` and `m ≤ 510`, so the true
value `w²·m` can reach `~2^71` — outside `UInt64`. `w * w` still fits
(`≤ 2^62`), so the significand square is a native `UInt64` multiply;
only the final scaling by `m` needs `widemul(::UInt64, ::UInt64) →
UInt128`. `widemul` compiles to a single native widening multiply on
x86-64 and AArch64 under Julia 1.12.

**Correction proof (accounts for integer flooring, `q = Q ÷ m`, not
`Q/m`).** Let `q_r = Q/m` be the exact rational quotient and
`q = ⌊q_r⌋ = Q ÷ m` the integer quotient. Two elementary facts:

1. `sqrt(q_r) − sqrt(q) ≤ sqrt(q_r) − sqrt(q_r − 1)`. For `q_r > 1/4`
   (equivalently `Q > m/4`, true here since `Q = 2^62` and `m ≤ 510`),
   `sqrt(q_r) − sqrt(q_r − 1) < 1/(2·sqrt(q_r − 1)) < 1`. Therefore
   `⌊sqrt(q_r)⌋ − ⌊sqrt(q)⌋ ∈ {0, 1}`.
2. `w = isqrt(q) = ⌊sqrt(q)⌋` satisfies `w² ≤ q ≤ q_r`, hence
   `w²·m ≤ Q` — the downward correction condition
   `w²·m > Q` is false, so the downward step never runs.
3. The true floor is `w* = ⌊sqrt(q_r)⌋ ∈ {w, w+1}` by fact 1. If
   `w* = w + 1` then `(w+1)² ≤ q_r`, i.e. `(w+1)²·m ≤ Q`, so the upward
   step fires exactly once. If `w* = w` then `(w+1)² > q_r`, i.e.
   `(w+1)²·m > Q`, so the upward step does not fire.

Hence a single conditional correction step in each direction suffices;
no `while` loop is needed, and the debug assertions above confirm the
correction bound at every call. Result has ≥ 27 significant bits.
Exhaustive T3 coverage per format enumerates both the fires-once and
fires-not cases, so the proof is checked, not trusted.

### 4.8 Hypot — exact sum, then sqrt of an **inexact** window

`x² + y²` is a two-term sum of exact 16-bit-sig products — `window_sum`
applies, and both terms are non-negative: no cancellation, premise safe.
But the square root then acts on an `Aligned` that may carry a tail, and
`√(w + ε)` can cross a boundary `√w` does not. The decidable algorithm:

1. `a = window_sum(x², y²)` — exact-or-tailed, non-negative.
2. Parity-shift `W = a.w << sh` (`sh ∈ {0,1}`; `a.w ≤ 2^57`, so
   `W ≤ 2^58`, safe), `s = isqrt(W)`, `r = W − s²`.
3. Enclose the true root: lower `s` (strict unless `r == 0 && !a.tail`),
   upper `s + 1`.
4. Project both endpoints under the caller's `ρ` with matching sidedness
   (§7.5). Equal result codes → done. This is the common case: the
   enclosure is ~2⁻²⁰-relative wide, rounding cells ~2⁻¹⁰.
5. If the endpoints straddle exactly one boundary `b` (a dyadic with
   ≤ 11-bit significand — candidate datum or midpoint), decide by an
   **exact integer comparison of the sum against `b²`**. `b²` has a
   ≤ 22-bit significand and an exponent `2·e_b` that can range up to
   `2·emax ≈ 300` binades from `a.w`'s units, so the aligned `b²` does
   not fit in a single `UInt64`. Represent aligned `b²` as an
   `NTuple{2,UInt64}` (low/high limbs) or a `UInt128`, and compare
   against `(a.w, 0)`/`UInt128(a.w)` respectively: differing integer
   parts differ by ≥ 1 unit, which the tail (< 1 unit) cannot flip;
   equal integer parts are decided by `a.tail` (`true` ⇒ strictly
   above `b`; `false` ⇒ exactly `b`, apply the mode's exact-point
   rule). Always decisive, always integer, never iterates.

This is §9's agree-or-refine loop collapsed to one guaranteed refinement,
allocation-free.

### 4.9 Convert

`unpack` source code point → `Datum` → hand to Projection for the target
format. Exact by construction; no window arithmetic involved.

### 4.10 Exports of this file

`Datum`, `Aligned`, `AlignedX`, `ExactResidue`, `AlignedResidue`,
`cmp_mag`, `window_sum`, `window_mul`, `window_fma`, `window_faa`,
`window_div`, `window_div_exact`, `window_sqrt`, `window_rsqrt`,
`window_hypot_sum` (steps 1–3 in §4.8; boundary resolution lives with
Projection, which owns sidedness). Nothing else. No format knowledge,
no policy knowledge, no `FloatByte`.

### 4.11 Which carrier flows where — the routing rules table

Every kernel selects the narrowest carrier its operation is proved to
admit; the projector reads whichever arrives. The routing is compile-
time once the projection's `R` (rounding-mode) type parameter is known,
so there is no runtime dispatch on carrier type in the hot path.

| operation | deterministic + `N ≤ 32` stochastic | `N > 32` stochastic |
|:---|:---|:---|
| `Add`, `Subtract` | `Aligned` | `AlignedX` |
| `Multiply` | `Datum` (exact) | `Datum` |
| `FMA` | `Aligned` (2-term after exact product) | `AlignedX` |
| `FAA` | `Aligned` — see §4.4 | `AlignedX` |
| `Sqrt` | `Aligned` (`isqrt` remainder is the tail bit) | `Aligned` |
| `RSqrt` | `Aligned` (integer enclosure, §4.7) | `Aligned` |
| `Hypot` | `Aligned` with §4.8 boundary refinement | `Aligned` |
| `Divide`, `Recip` | `Aligned` (49-bit quotient decides `N ≤ 32`) | **`ExactResidue`** |
| reduction handoff (§13.1) | `AlignedResidue{Lhi}` | `AlignedResidue{Lhi}` |

Projection carries four `project_code` methods (§7.1):

```julia
project_code(::Type{FR}, ρ, x::Aligned,             entropy)::UInt8
project_code(::Type{FR}, ρ, x::AlignedX,            entropy)::UInt8
project_code(::Type{FR}, ρ, x::ExactResidue,        entropy)::UInt8
project_code(::Type{FR}, ρ, x::AlignedResidue{Lhi}, entropy)::UInt8 where {Lhi}
```

The `Aligned` method is the hot pointwise path. The others are called
only when the kernel proved the extra precision is needed (`AlignedX`
for stochastic `N > 32` and FAA), the residue is non-dyadic
(`ExactResidue` for Divide with `N > 32`), or a reduction closes
(`AlignedResidue` from `LimbAcc`/`MulAcc`). Stage 4 verifies each
method against `reference/` on its declared domain, and the
non-`Aligned` methods are additionally verified against `Aligned` on the
common domain where both apply. Any carrier whose bound is not proved by
Stage 4 routes to the cold arbitrary-precision adapter in `reference/`;
no production carrier silently widens.

---

## 5. The Byte Domain — `domain.jl`, `formats.jl`, `show.jl`

### 5.1 The type

```julia
abstract type Signedness end
struct SignedFormat   <: Signedness end
struct UnsignedFormat <: Signedness end
abstract type DomainKind end
struct FiniteDomain   <: DomainKind end
struct ExtendedDomain <: DomainKind end

primitive type FloatByte{K,P,S<:Signedness,D<:DomainKind} <: AbstractFloat 8 end
```

The marker names are deliberately long: `Signed` and `Unsigned` are Base
abstract integer types this package also dispatches on; rebinding them
inside the module is a defect factory. Users never type the markers — the
120 exported aliases carry the draft's names:

```julia
# formats.jl — 120 const aliases following draft §4.2 naming, written out
# by a plain loop in the source generator (not eval at load):
const Binary8p4se = FloatByte{8,4,SignedFormat,ExtendedDomain}
# ... all 120, plus exports.
```

One byte, concrete, `isbits`. Verify in Stage 2 (T1): `isbitstype`,
`sizeof == 1`, `Base.elsize(Vector{T}) == 1`, `reinterpret(UInt8,
::Vector{T})` zero-copy. Because the type is concrete, "format types are
propagated, never rebuilt" is vacuous — the defect class does not exist.

`validate_format(F)` checks `K ∈ 3:8`, the `P` range for `Σ`, at every
public seam that first receives a format type; it is written so constant
formats fold the check away entirely (G9 covers this). No mutable format
registry exists.

Subtyping `AbstractFloat` is retained (ComposedPlan decision; the plan
verified layout and zero-copy properties). The phase-led guide's caution —
Base's `AbstractFloat` generic fallbacks carry promotion/rounding
expectations — is honored by **staging**: until Stage 5 lands the deliberate
surface, the package defines only what it has tested, and G10 (surface
totality) plus explicit `MethodError`s keep untested generic fallbacks
from silently producing host-float semantics.

### 5.2 Traits and format queries — all fold at compile time

Internal traits, one-line formulas from §1:

```julia
kbits(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = K
pbits(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = P
issigned(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = S === SignedFormat
isextended(::Type{FloatByte{K,P,S,D}}) where {K,P,S,D} = D === ExtendedDomain
bias(T) = issigned(T) ? 1 << (kbits(T)-pbits(T)-1) : 1 << (kbits(T)-pbits(T))
# expbits, trailingbits, code_nan, code_pinf, code_ninf, code_one,
# maxfinite_code, minpositive_code, emax — same pattern; absent specials
# return nothing.
```

Public queries follow Julia style with a documented mapping to the draft's
names ([§5.14](IEEE_D1.md#L1975)); all accept `::Type{F}` (never require
an instance):

| exported | draft |
|:---|:---|
| `bitwidth(F)` | `BitwidthOf` |
| `Base.precision(F)` | `PrecisionOf` |
| `signedness(F)`, `domain(F)` | `SignednessOf`, `DomainOf` |
| `exponent_bitwidth(F)`, `trailing_significand_bits(F)` | `ExponentBitwidthOf`, `TrailingSignificandBitwidthOf` |
| `exponent_bias(F)` | `ExponentBiasOf` |
| `maxfinite(F)`, `minfinite(F)`, `minpositive(F)`, `maxsubnormal(F)`, `minnormal(F)` | `MaxFiniteOf`, … (value queries return the datum) |

Gate **G9**: every trait and query constant-folds at every alias
(`@inferred` plus code inspection asserting literal constants).

### 5.3 Code identity

```julia
fromcode(::Type{T}, code::Integer) where {T<:FloatByte} =
    0 <= code < (1 << kbits(T)) ? reinterpret(T, UInt8(code)) :
        throw(ArgumentError("code $code out of range for $(T)"))
codepoint(x::FloatByte) = reinterpret(UInt8, x)
isvalidcode(::Type{T}, code::Integer) = 0 <= code < (1 << kbits(T))
```

`fromcode` is the only safe code-identity input. `T(value::Real)` means
numeric value at every integer width — there is no `Unsigned`-means-code
rule, and the raw and mathematical constructors can never be confused:
construction from a real value is visibly projecting
(`quantize(T, value, ρ; entropy)`; `T(value)` uses the fixed package
default and documents it).

`reinterpret` into a `FloatByte` appears in exactly two places in `src/`:
these definitions, and the byte-view lowering in `execution/dense.jl`. A
test greps for it (§17.6). It is documented unsafe for users: it can
fabricate nonzero high bits when `K < 8`, and every kernel relies on the
canonical-value invariant.

### 5.4 The decode atlas

One `const` `Memory{Float64}` of 13 296 entries in `format_id` order
(106 368 B), plus the `Float32` twin (53 184 B) — `Memory` rather than
`Vector`: fixed-size, one less indirection, and it documents immutability
of intent. Built by a plain generator at precompile time, checked against
`reference/` in T1:

```julia
decode(x::FloatByte{K,P,S,D}) where {K,P,S,D} =
    DECODE64[FORMAT_OFFSETS[format_id(K,P,S,D)] + Int(codepoint(x)) + 1]
```

One const-`Memory` load for `FORMAT_OFFSETS`, one for `DECODE64`, two adds
and a bounds-elided indexed read; `format_id` folds so the offset index is
a constant per alias. The Stage 2 benchmark compares the atlas against
per-format 256-entry tables (lookup latency, image size, load time); the
loser is deleted.

`unpack(x)::Datum` is the internal exact form: specials excluded by the
caller, then §1.3's field extraction into `(sig, exp, neg)` with
`sig = S`, `exp = E + 1 − P`. Pure bit arithmetic — it must not read the
atlas.

### 5.5 Ordering

Total order per [§5.12.1](IEEE_D1.md#L1851): **NaN first, below `−Inf`** —
the opposite end from Float64 sorting. Implementation: per-format
256-entry `UInt16` order-key table (`Memory{UInt16}`), key 0 reserved for
NaN, built once from the sorted decoded lattice. `TotalOrder`, `isless`,
sorting, and the min/max family consult keys, never floats. Counting
sort over keys above a measured length threshold
`SORT_COUNT_THRESHOLD` (a `const Int`, set at Stage 6 from the
benchmark harness's measured crossover on the target microarchitecture
and pinned in `benchmark/thresholds.md`); `Base.sort` below the
threshold; caller-owned `SortWorkspace` for repeated zero-allocation
sorts — no hidden shared scratch.

### 5.6 `show` discipline

`show` never throws and never projects: format alias, code point, decoded
value from the atlas (`Binary8p4se(0x42 ↦ 1.5)` style; exact decimal since
every datum is a Float64-representable dyadic). A failing testset must be
able to print anything, including invalid bytes fabricated by
`reinterpret` — `show` masks and flags out-of-range codes rather than
erroring. Docstrings: no `\(`, `\[`, `\text` — invalid string escapes
that can render a package unloadable; the §17.6 grep-enforced structural
rules include a docstring-escape row that fails CI on any occurrence in
`src/`.

---

## 6. `policies.jl`

```julia
abstract type RoundingMode end
struct NearestTiesToEven <: RoundingMode end
struct NearestTiesToAway <: RoundingMode end
struct TowardPositive    <: RoundingMode end
struct TowardNegative    <: RoundingMode end
struct TowardZero        <: RoundingMode end
struct ToOdd             <: RoundingMode end
struct StochasticA{N}    <: RoundingMode end     # N random bits, §5.7.4
struct StochasticB{N}    <: RoundingMode end
struct StochasticC{N}    <: RoundingMode end

abstract type SaturationMode end
struct SatFinite    <: SaturationMode end
struct SatPropagate <: SaturationMode end
struct SatNone      <: SaturationMode end

struct Projection{R<:RoundingMode, S<:SaturationMode} end   # ρ
const DEFAULT_PROJECTION = Projection{NearestTiesToEven, SatNone}()

# Entropy class of a rounding mode. Determines whether a kernel consumes
# the entropy word and, downstream, which plan classes and SIMD forms are
# admissible (§11.4, §12.1). A trait, not a type hierarchy: every
# RoundingMode subtype gets exactly one method here, checked exhaustively
# by G10.
abstract type EntropyClass end
struct Deterministic <: EntropyClass end
struct Stochastic    <: EntropyClass end

entropy_class(::Type{<:RoundingMode})  = Deterministic()
entropy_class(::Type{<:StochasticA})   = Stochastic()
entropy_class(::Type{<:StochasticB})   = Stochastic()
entropy_class(::Type{<:StochasticC})   = Stochastic()

# Signature-level convenience: forwards to the R parameter.
entropy_class(::Type{Signature{Op,Fx,Fy,Fr,R}}) where {Op,Fx,Fy,Fr,R} =
    entropy_class(R)
```

Policy is dispatch, never a branch (ComposedPlan invariant 16).
`RoundOf`/`SatOf` ([§5.15](IEEE_D1.md#L2005)) are accessors on
`Projection`.

Stochastic modes are **stateless singletons parameterized by `N` only**;
the random word travels separately (§12). Storing bits inside the mode
value would conflate policy identity with per-call entropy and break the
one-word-per-logical-output contract. Validation kept from the phase-led
guide: `N > 0` checked at mode construction; the entropy word checked
`< 2^N` at the seam. Kernels never touch global RNG state — results are
reproducible and exhaustively testable by construction.

The conformance-oriented interface **requires** the explicit `ρ` argument;
only documented convenience wrappers use `DEFAULT_PROJECTION`.

---

## 7. `projection.jl` — the only code-producing implementation

```julia
project_code(::Type{Tr}, ρ, a::Aligned, entropy::UInt64)::UInt8
project_interval(::Type{Tr}, ρ, lo::Aligned, hi::Aligned, entropy::UInt64)::UInt8
```

Internally exactly `round_to_precision`, then `saturate`, then `encode` —
never merged: rounding precedes saturation, the rounding mode remains
visible to saturation, and saturation never rounds
([§4.7.3](IEEE_D1.md#L665) NOTE 1, [§4.7.5](IEEE_D1.md#L770); same
draft-line anchors as invariant 2 in §2). Projection knows formats and
policies; it knows nothing of operations, arrays, tables, or plans.

### 7.1 Round to precision — unbounded exponent first

Rounding to `P` significant bits happens with **unbounded** exponent range
(subnormal quantum clamp applied, overflow left in place); saturation then
classifies. This ordering makes the NTE half-ulp-above-MaxFinite rule
(concepts §15, law 5) fall out with no special case.

```julia
function round_to_precision(a::Aligned, P, B, mode, entropy)  # -> (sig, E, neg) | :zero
    h      = 63 - leading_zeros(a.w)                # w == 0 handled before this
    msbexp = a.exp + h
    ulpexp = max(msbexp, 1 - B) + 1 - P             # subnormal quantum clamp
    k      = ulpexp - a.exp                         # ulp bit position within w
    # k ≥ 33 always: msbexp ≥ a.exp + 40 (window guarantee) and P ≤ 8.
    keep   = k >= 64 ? UInt64(0) : a.w >> k
    rbit, sticky = split_fraction(a.w, k, a.tail)   # bit k-1; OR of bits below + tail
    keep  += round_away(mode, keep, rbit, sticky, a.neg, entropy, P) ? 1 : 0
    if keep == UInt64(1) << P                        # carry out of the precision
        keep >>= 1; ulpexp += 1
    end
    keep == 0 ? :zero : (keep, ulpexp + P - 1, a.neg)   # canonical (S, E)
end
```

The `k ≥ 33` bound is the window guarantee doing its job: the fraction
below the ulp is known to ≥ 33 bits plus a sticky, which decides every
deterministic mode and every stochastic mode with `N ≤ 32`.

**`project_code` for the exact-residue carriers.** At
`R = Stochastic*{N > 32}`, kernels emit `AlignedX` (dyadic residue) or
`ExactResidue` (non-dyadic Divide/Recip residue) and the projection
dispatches to carrier-specific methods:

```julia
# Dyadic residue: the exact dropped remainder already lives in the frame
# rem.sig at rem.exp. Compare the top-N bits of the residue-in-window
# against the entropy word as pure integer arithmetic.
function project_code(::Type{Tr}, ρ::Projection{R,S},
                      a::AlignedX, entropy::UInt64)::UInt8 where {Tr, R<:Union{StochasticA,StochasticB,StochasticC}, S}
    # ... use a.a for the coarse decision; a.rem furnishes the top-N tail bits;
    # threshold  = entropy >> (64 - N); compare, then encode via saturate∘encode.
    # N ∈ 1:63 by construction: _stochastic_bits(R) is a compile-time
    # constant and §16 places N == 64 out of the covered stochastic range.
end

# Non-dyadic residue: floor(widemul(num, 2^N) / den) is the top-N-bit
# fraction the predicate needs, computed with a UInt128 dividend.
function project_code(::Type{Tr}, ρ::Projection{R,S},
                      r::ExactResidue, entropy::UInt64)::UInt8 where {Tr, R<:Union{StochasticA,StochasticB,StochasticC}, S}
    N          = _stochastic_bits(R)                   # compile-time constant
    @assert 1 <= N <= 63                                # debug: shift-by-N safety
                                                        # (N == 64 is not a defined mode; see §16)
    scaled     = widemul(r.num, UInt64(1) << N)         # UInt128; UInt64(1) << N well-defined for N ≤ 63
    top_n      = UInt64(div(scaled, UInt128(r.den)))    # 0 .. 2^N - 1
    threshold  = entropy >> (64 - N)                    # top-N bits of the entropy word
    up         = top_n >= threshold                     # StochasticA form; B/C differ per §5.7.4
    # ... combine with r.q, r.qexp, r.neg; then saturate ∘ encode.
end
```

Deterministic and `N ≤ 32` stochastic never see these carriers — they
use `window_div` → `Aligned` and the existing `project_code(::Aligned)`
method above. `project_code(::AlignedResidue{Lhi})` (§13.1's reduction
handoff) reads `hi[Lhi]`'s leading nonzero as the rounding decision's
coarse magnitude, `hi[1..Lhi-1]` for the below-boundary bits that place
the round bit, and `lo` as the sticky OR of every lower limb; the exact
formula follows the same shape as `Aligned` (round bit + sticky) but
generalized over the multi-limb integer, and `AlignedResidue` never
appears in the pointwise fast path.

**Rounding predicate reference — uniform expression.** With
`ν = discarded_fraction / discarded_scale` the exact rational below the
ulp (extracted per carrier: bit slice from `Aligned`, `Aligned + rem` for
`AlignedX`, `num/den` for `ExactResidue`, `hi[1..Lhi-1] + lo`-sticky for
`AlignedResidue`), and `M = 2^N`, `0 ≤ R < M` the entropy word's top-`N`
bits:

| mode | round away from zero when |
|:---|:---|
| `TowardZero` | never |
| `TowardPositive` | `ν > 0` and `!neg` |
| `TowardNegative` | `ν > 0` and `neg` |
| `NearestTiesToAway` | `2ν ≥ 1` (equivalently `2r ≥ d` in integer form) |
| `NearestTiesToEven` | `2ν > 1`, or tie and retained code is odd (**code-point evenness** per §5.7.4's `CodeIsEven` for `P == 1`) |
| `ToOdd` | `ν > 0` and retained code is even (bit-OR the kept LSB, never increment) |
| `StochasticA{N}` | `floor(ν·M) + R ≥ M` |
| `StochasticB{N}` | `floor(ν·2M) + (2R + 1) ≥ 2M` |
| `StochasticC{N}` | `RNITE(ν·M) + R ≥ M` |

Every predicate evaluates as pure integer arithmetic. Never convert `ν`
to `Float64`, never build a `BigFloat`, and never introduce a
computation whose result depends on host rounding mode. `N` is validated
at mode construction (`1 ≤ N ≤ 63`), and `R < 2^N` is asserted at the
seam — both are compile-time checks under a typed `IndexedEntropy`.

`round_away(mode, …)` is the internal per-singleton dispatch that
realizes the predicate table above; one method per rounding-mode
singleton, no runtime branches:

| mode | away iff |
|:---|:---|
| `NearestTiesToEven` | `rbit && (sticky \|\| isodd(keep))` — with `P == 1` evenness defined on the **code point** per the draft's `CodeIsEven` ([§5.7.4](IEEE_D1.md#L704)): parity depends on zero and the parity of `Q + B`, not the integer significand. Make it a named helper, test it independently, pin it at T2 — do not guess parity semantics |
| `NearestTiesToAway` | `rbit` |
| `TowardPositive` | `!neg && (rbit \|\| sticky)` |
| `TowardNegative` | `neg && (rbit \|\| sticky)` |
| `TowardZero` | never |
| `ToOdd` | truncate, then force the kept LSB odd when inexact — a **bit-or**, not an increment; implement the draft's RoundAway table literally |
| `StochasticA/B/C{N}` | compare the top-`N` fraction bits (with sticky) against the entropy word per the draft's three predicates (`StochasticC` via the `RNITE` helper), as pure integer comparisons — never convert the fraction to Float64 first. `N ≤ 32`: decided exactly from the window. **`N > 32`**: the window's ~33-bit fraction knowledge is insufficient; kernels route the exact residue to a carrier-specific `project_code` method — `AlignedX` for Add/Multiply/FMA/FAA/Sqrt/RSqrt (dyadic residue), `ExactResidue` for Divide/Recip (non-dyadic residue with `den = y.sig`), and the enclosure-loop path (§9.2) for Hypot. Each method compares against the entropy word as pure integer arithmetic: `AlignedX` uses the exact dyadic tail directly; `ExactResidue` computes `floor(widemul(num, UInt64(1) << N) / den)` (a `UInt128` value) and compares against `entropy << (64 - N)`. `N` is bounded to `1:63` at mode construction; `N == 64` would need a two-word widemul and is not a requested mode. T5 exhausts all `R` for small `N` (including fractions just below/at/above one half) and pins maximum-`N` edges |

### 7.2 The reference rounding algorithm (differs by construction)

`reference/rational_projection.jl` decides the same questions from exact
`Rational{BigInt}` quotient/remainder — floor `q` at quantum
`2^Q, Q = max(e, 1−B) − P + 1`, then:

| mode | increment iff |
|:---|:---|
| toward zero | never |
| toward positive / negative | remainder ≠ 0 and sign matches |
| nearest ties away | `2r ≥ d` |
| nearest ties even | `2r > d`, or `2r == d` and `q` odd (`P == 1` via the code-parity helper) |
| to odd | remainder ≠ 0 and `q` even |

The production path uses window guard/round/sticky bits; the reference
uses exact quotient/remainder. **They must not share the rounding decision
code** — their exhaustive agreement is meaningful precisely because the
algorithms differ (oracle discipline, §15). Shared format traits and
special-code constants are fine.

### 7.3 Saturate

Input: rounded `(sig, E, neg)` with unbounded `E`; target traits; the
rounding mode. An ordered decision table mirroring the draft's first-match
behavior; magnitude compared against `MaxFiniteOf(Tr)` as `(sig, E)`
pairs, never as floats:

- in range → pass through;
- `SatFinite` → clamp to `M_hi`/`M_lo`;
- `SatPropagate` → clamp finite overflow; preserve representable
  infinities (infinite **inputs** never reach here — code-point rows in the
  kernels consult the saturation mode);
- `SatNone` → per [§4.7.5](IEEE_D1.md#L770): overflow to `±Inf` when
  `Δ = Extended` and the mode rounds away; to `M_hi`/`M_lo` under a
  directed mode pointing inward; NaN where the draft says NaN (finite
  domain, negative in unsigned, …). Implement the draft's table
  literally.

T2 enumerates the full cross-product: 3 saturation modes × every rounding
mode × signed/unsigned × finite/extended × {NaN, ±Inf, exact extrema, one
step beyond, a finite value that rounds beyond}. The four cases the
phase-led guide singled out for protection — directed rounding forcing an
extremal finite under `SatNone`; negative overflow in unsigned extended
becoming NaN; infinities propagating only where representable; NTE
rounding just-above-`M_hi` back to `M_hi` — are named test cases, not
prose.

### 7.4 Zero and encode

`w == 0 && !tail` → code 0 (one zero, no `−0` bookkeeping anywhere;
external `−0.0` collapses at the bridge, §14.3). `w == 0 && tail` cannot
reach projection from any tier-1 kernel (§4.2); assert it away. Encode is
§1.4, reached only from saturate with in-range canonical arguments.

### 7.5 Sidedness

`project_interval` projects the lower endpoint with **positive** sidedness
and the upper with **negative**, same `ρ`, same entropy word. Agreement ⇒
that code; disagreement ⇒ caller-specific resolution (exact boundary
comparison for Hypot §4.8; precision raise for transcendentals §9;
`InconclusiveProjectionError` at the ceiling).

---

## 8. Operation kernels — `operations/`

### 8.1 The kernel shape

```julia
function op_kernel(sig::Signature{AddOp,Fx,Fy,Fr,R}, e::UInt64, cx::UInt8, cy::UInt8)::UInt8
    # 1. SPECIAL ROWS — pure code-point tests against the draft's pattern
    #    tables (§5.3.2, "first match wins" is normative), each row a cited
    #    draft line. NaN in → NaN out; Inf rows per op; zero rows per op.
    # 2. unpack survivors (finite, nonzero where the op requires)
    # 3. one window routine from §4
    # 4. one projection: project_code(Fr, ρ, aligned, e)
end
```

The special rows **are** the spec: the draft states them as code-point
patterns, and code-point tests are one indirection closer to the text than
decoding and asking `isnan`. Keep the draft's row order. Do not assume
IEEE-754 host behavior anywhere — the differences are encoded explicitly
and each is a named test (§8.2 notes column). Table-driven tests are
generated from each operation's ordered case list, then the finite kernel
is tested separately, so an omitted row is visible rather than silently
delegated to host math.

### 8.2 Per-operation notes (algebraic)

| op | draft | kernel core | notes |
|:---|:---|:---|:---|
| `Add`/`Subtract` | [§5.10.3](IEEE_D1.md#L1008) | `window_sum` | `Inf + (−Inf) = NaN` row before unpack |
| `Multiply` | [§5.10.4](IEEE_D1.md#L1052) | `window_mul` | `0 × Inf = NaN` row |
| `Divide` | [§5.10.5](IEEE_D1.md#L1093) | `window_div` | **`x/0 = NaN`** (not Inf), `x/±Inf = 0` — draft NOTEs 1–2; easy to get IEEE-754-wrong |
| `FMA` | [§5.10.6](IEEE_D1.md#L1131) | `window_fma` | specials via the equivalent `Add(Multiply(x,y),z)` rows; single projection |
| `FAA` | [§5.10.7](IEEE_D1.md#L1182) | `window_faa` | three-term §4.4; associativity property test; cancellation + large exponent gaps mandatory in T3 |
| `Sqrt` | [§5.10.8](IEEE_D1.md#L1232) | `window_sqrt` | negative → NaN; `Sqrt(+Inf) = +Inf` |
| `Recip` | [§5.10.8](IEEE_D1.md#L1232) | `window_div` | `Recip(0) = NaN`, `Recip(±Inf) = 0` |
| `RSqrt` | [§5.10.8](IEEE_D1.md#L1232) | `window_rsqrt` | §4.7 direct form |
| `Hypot` | [§5.10.14](IEEE_D1.md#L1492) | §4.8 | algebraic, **not** a transcendental; Inf rows dominate NaN rows here — read the pattern table, don't assume |
| `Abs`/`Negate`/`CopySign` | [§5.10.1](IEEE_D1.md#L945)/[.2](IEEE_D1.md#L976) | code-point bit ops | sign-bit flip/copy on signed; unsigned negation follows the draft's rows (NaN unless zero) |
| `Convert` | [§5.9.1](IEEE_D1.md#L918) | unpack → project | the only cross-format tier-1 op |

### 8.3 Selection family

`Minimum/Maximum` + eight variants ([§5.11.1–.3](IEEE_D1.md#L1602)),
`Clamp` ([§5.11.4](IEEE_D1.md#L1737)): order-key comparisons plus each
variant's NaN/infinity quieting rules. Same-format selection returns an
operand code unchanged — no projection (the result is an input datum;
re-projecting would violate "one projection"). Cross-format variants go
through Convert semantics where the draft says so.

### 8.4 Comparisons, predicates, class

[§5.12](IEEE_D1.md#L1778)–[§5.13](IEEE_D1.md#L1895): Booleans/enums, stop
after decode — implemented as order-key/code tests, no projection, no
entropy. `Class` partitions every code point into exactly one of the eight
classes; classification is computed from code fields and *proved
equivalent to decoded classification exhaustively* (T1). These ops bypass
`defined_code` (they produce no code) through their own seam
`defined_predicate`; both seams are catalogue-generated and G10-covered.

### 8.5 Neighbours

`NextGreaterThan`/`NextLessThan` ([§5.16](IEEE_D1.md#L2027)): order-key
walks with explicit terminal rows (`NextGreaterThan(+Inf) = NaN`, unlike
IEEE-754 `nextUp`). T1 tests the whole graph per format: every nonterminal
step decodes strictly monotone; the two functions are mutual inverses away
from terminals; signed transitions cross the unique zero; nothing wraps at
`0` or `2^K − 1`.

---

## 9. Transcendentals — `operations/elementary.jl`

25 unary operations plus `ArcTan2`/`ArcTan2Pi`
([§5.10.9–.13, .15–.16](IEEE_D1.md#L1273)). Hypot is **not** here.

### 9.1 Always-on unary tables

Every unary transcendental at every format ships as a 256-byte table per
requested `(op, format, ρ)` **for deterministic ρ only** — built on
first use, then immutable (`Memory{UInt8}`). A table entry **is** a call
to `defined_code`; tables never define behavior. A table becomes visible
only complete and checksummed (SHA-256 recorded beside it).

Stochastic `ρ` on a transcendental never tabulates: a table entry would
freeze one entropy decision per code point, which contradicts the
per-invocation entropy contract (§11.2, §12.1). Stochastic
transcendentals route to the §9.2 enclosure loop per call. G-gate: the
table-build path grep-rejects any `Projection{<:Stochastic*, ...}` at
plan preparation.

### 9.2 The rigorous enclosure loop (cold path)

A fixed-precision `BigFloat` evaluation is not a correct-rounding
argument. The loop is Ziv-style with a certified directed enclosure:

```text
1. exact/special cases symbolically (draft rows; exact points like
   Log(1) = 0 short-circuit — the per-op exact-case table is finite and cited)
2. directed MPFR enclosure at working precision p (start p = 32):
   lo ≤ f(x) ≤ hi, endpoints dyadic, sidedness flags carried
3. project lo with positive sidedness, hi with negative, same ρ, same
   entropy word (§7.5)
4. equal codes -> return; else p *= 2, repeat
5. at the caller-supplied ceiling: throw InconclusiveProjectionError
   carrying (signature, operand code, precisions, final enclosure) —
   never a midpoint, never a guess
```

Asymptote sidedness is symbolic: `Tanh(large)` is `1 − ε`, never `1`.
If an enclosure cannot be certified for some function, that function is
**registered as approximate under a distinct name** (§16) — it does not
ship as a defined operation.

With ≤ 256 inputs per `(op, format)`, hard cases are findable: Stage 6
enumerates the precision each entry needed and records the per-`(op,
format)` maximum in a shipped `Memory{UInt16}` constant table
(`P_CEILING_TABLE :: Memory{UInt16}`, indexed by
`(op_id, format_id)`), so production builds start at known-sufficient
precision and the loop is one pass in practice. The published table is
part of the release manifest and covered by its own SHA-256 digest.

Stage 6's exit rule refuses to ship any `(op, format)` whose measured
ceiling meets or exceeds `1 024` bits without an accompanying §16
approximate registration: at that precision the elementary function is
known-hard on that argument, and shipping it as exact would silently
switch users onto the cold Ziv-doubling path forever. The gate is
mechanical (read Stage-6 CSV, filter for `p >= 1024`, cross-check
`conformance_report()`); any unregistered entry blocks the release.

`ArcTan2` (≤ 65 536 inputs per format pair) uses the same loop; its table
policy follows the binary rule (§11.2). `ArcTan2(0,0) = NaN` and the
quadrant rows come first, as code-point rows.

MPFR appears in this file and in `reference/` — nowhere else (grep-
enforced, §17.6).

---

## 10. The seam — `semantics.jl`, `catalogue.jl`

### 10.1 Operations, signatures, `defined_code`

```julia
abstract type Operation end
struct AddOp <: Operation end;  struct MultiplyOp <: Operation end
struct FMAOp <: Operation end   # ... one singleton per catalogue op

struct Signature{Op, Fx, Fy, Fr, R}  end     # all type parameters; isbits
# arity variants: {Op,Fx,Fr,R} unary; {Op,Fx,Fy,Fz,Fr,R} ternary

defined_code(sig::Signature, entropy::UInt64, codes::UInt8...)::UInt8
```

`defined_code` is the internal seam: scalar wrappers, table builders,
dense loops, block lanes, κ measurement, and conformance all call this
function and nothing lower. A table entry is a call to it, so no
correctness reasoning remains at any use site.

The public seam wraps it:

```julia
apply(op::Operation, ::Type{Fr}, ρ::Projection, xs...)   # heterogeneous operands
add(Fr, ρ, x, y)  = apply(AddOp(), Fr, ρ, x, y)          # thin discovery wrappers
fma(Fr, ρ, x, y, z) = apply(FMAOp(), Fr, ρ, x, y, z)
```

`apply` centralizes arity checks and format resolution; operands supply
`formatof`/`decode` methods (external-format values enter here, §14.3);
an unsupported operand is a `MethodError`, never a silent fall-through to
host arithmetic. No public unary/binary/ternary class hierarchy — arity
is an implementation trait.

Entropy contract: deterministic modes ignore the word; stochastic modes
consume exactly the word passed, exactly once, including for exact and
special results — fixed consumption is what makes slicing and replay
auditable. `entropy = 0` under a stochastic mode is legal input, not a
sentinel.

### 10.2 `catalogue.jl`

One `const` vector of a plain `OpSpec` struct — metadata only: op name,
arity, operand domains, whether it projects, draft section, conformance
group (Annex G), table policy class. Generated from it: exported
wrappers, `apply` methods, G10 surface tests, conformance declaration
strings, release-tier loops. Semantic bodies are ordinary reviewed Julia
in `operations/`; the catalogue never contains behavior.

---

## 11. Plans, tables, execution

### 11.1 Three plan classes (Stage 6)

```julia
DenseComputePlan{Sig}     # tier-1 kernel over a byte view
UnaryTablePlan{Sig}       # the always-on 256-byte table
GenericArrayPlan{Sig}     # stride/view/odd-axis correctness adapter
```

Immutable; own their tables strongly (the plan holds the `Memory{UInt8}`
reference — eviction cannot invalidate a running plan). The planner is
cost-only: plan choice may never change a byte. `plan_traits(plan)`
reports workspace needs and fallbacks **before** execution. Later classes
(`BinaryTablePlan`, SIMD, threaded, direct-code laws) arrive one at a
time in Stage 11, each with a measured win and an exhaustive equivalence
certificate, else deleted.

### 11.2 Table policy (ComposedPlan §10.2, adopted)

| candidate | policy |
|:---|:---|
| unary, any `K` (≤ 256 B) | always — L1-resident |
| binary, `K ≤ 5` (≤ 1 KiB) | yes — L1-resident |
| binary, `K ∈ 6:8` (4–64 KiB) | only on measured reuse; break-even ≈ 3× the array ≈ 200 k elements per signature |
| ternary, any `K` | never — the kernel beats what the build amortizes |
| stochastic ρ as a result table | forbidden — a table entry would freeze an entropy decision. Applies to §9 transcendental tables as well: transcendentals ship 256-byte tables only for deterministic ρ; stochastic ρ for a transcendental takes the enclosure-loop cold path (§9.2) per call |
| artifact packs | deferred until a demonstrated build-cost problem |

Binary indexing: `table[(Int(ca) << Ky) | Int(cb) + 1]`, row-major in the
same `format_id` conventions as everything else.

### 11.3 Build discipline

allocate → fill every entry via `defined_code` → SHA-256 digest → publish.
Never visible partially built. Generation metadata (package version,
draft revision, signature) recorded with the digest. Build-cost target:
**format-independent** — ≤ 80 µs / ≤ 44 allocations for a 64 K binary
table on 8 threads, wide-spread formats included, because no escalation
path exists to differ by.

### 11.4 Dense lowering

```julia
function execute!(plan::DenseComputePlan{Sig}, out, a, b)
    ca = reinterpret(UInt8, a); cb = reinterpret(UInt8, b)
    co = reinterpret(UInt8, out)
    @inbounds for i in eachindex(co, ca, cb)
        co[i] = defined_code(Sig(), entropy_for(plan, i), ca[i], cb[i])
    end
end
```

Byte views taken **once, outside the loop**; per-element there is no
registry lookup, format check, dictionary, lock, or strategy branch.
Aliasing `out` with an input is legal for same-format elementwise ops;
the generic plan detects overlapping strided views and copies, throwing
`AliasingError` only where semantics would differ.

No `@simd` on the generic form. Under a stochastic `R`, the body calls
`entropy_for` → `entropy_word` → four `mix64` invocations, and `@simd`
degrades to a no-op on any body with a non-inlineable state-carrying
call. Under a deterministic `R` the body is pure integer arithmetic and
vectorizable, but SIMD is a Stage 11 deliverable with its own equivalence
and feature-detection gates (§11.5); putting `@simd` here — where the
same plan class serves both `R` cases — either misleads the reader or
silently overpromises. Stage 11's SIMD plan is a separate plan class `DenseSIMDPlan{Sig}` whose
`execute!` method is defined only under a `where` clause that pins the
signature's entropy class to deterministic:

```julia
struct DenseSIMDPlan{Sig} <: AbstractPlan
    # ... table refs, tile size, feature flags ...
end

# Only defined for deterministic signatures. The `where`-clause enforces
# entropy class at method-resolution time: a stochastic-Sig call site is
# a genuine MethodError (produced by Julia's own dispatch), and the
# planner catches it and downgrades to DenseComputePlan. Do not manually
# construct a `MethodError`; a hand-thrown one carries the wrong CFG for
# `Base.showerror` and confuses tooling.
function execute!(plan::DenseSIMDPlan{Sig}, out, a, b
                  ) where {Sig <: OpSignature{Op,F} where {Op,F}}
    _ = entropy_class(Sig)::Deterministic   # compile-time constraint
    # @simd ivdep body ...
end
```

`entropy_class(Sig)` folds at compile time (G9), so the check compiles
away on the deterministic path and becomes a link-time exclusion on the
stochastic path. `@simd ivdep` is admitted in the body only after
inlining is proved by `code_llvm` inspection and the equivalence gate
passes (T4).

### 11.5 Threading and SIMD (Stage 11, gated)

Before threading anything, verify the serial path specializes —
SmallFloats' threading win was a missing function barrier in disguise (a
runtime `Symbol` → `Val` inside a builder). The benchmark preflight
aborts if a warm scalar path allocates. Threaded deterministic execution
must produce identical bytes. SIMD plans need exhaustive code
equivalence, feature detection at preparation, tail handling tested at
every short length, and a measured win; AVX2 has no byte gather usable
for a 64 KiB table.

---

## 12. Stochastic execution — `execution/stochastic.jl`

### 12.1 The indexed entropy word

```julia
struct IndexedEntropy
    seed::UInt64; stream::UInt64; invocation::UInt64
end

const ENTROPY_ALGO_VERSION = 0x01

# splitmix64 finalizer; canonical constants.
@inline function mix64(z::UInt64)
    z = xor(z, z >> 30) * 0xbf58476d1ce4e5b9
    z = xor(z, z >> 27) * 0x94d049bb133111eb
    xor(z, z >> 31)
end

# N is carried as a Val so `64 - N` folds; the shift is a constant-shift
# instruction, not a runtime-length shift. Kernels reach entropy_word via
# `entropy_word(ie, i, slot, Val(N))` where N comes from the signature's
# StochasticA/B/C{N} type parameter.
@inline function entropy_word(ie::IndexedEntropy, index::Int, slot::UInt32, ::Val{N})::UInt64 where {N}
    h = mix64(ie.seed + 0x9e3779b97f4a7c15 * ENTROPY_ALGO_VERSION)
    h = mix64(xor(h, ie.stream))
    h = mix64(xor(h, ie.invocation))
    h = mix64(xor(h, (UInt64(index) << 32) | UInt64(slot)))
    h >> (64 - N)                       # top N bits; (64 - N) is a compile-time constant
end
```

`@inline` on both is not decorative, and neither is the `Val{N}`: the
dense loop (§11.4) contains `defined_code` → kernel → `entropy_word` →
four `mix64` calls per element. Without inlining, the per-element cost
is dominated by the call chain. Without `Val{N}`, `h >> (64 - N)` is a
runtime-length shift (roughly 4× the latency of a constant shift on
x86-64) and the compiler cannot lift the shift out of loops over a
fixed-signature array. Stage 8's W3 check asserts zero call sites and
no runtime-length shifts remain in the lowered form of the stochastic
kernel.

Exact chaining frozen at Stage 8 under golden vectors; the requirement is
the shape — one mix per key component, top-`N` extraction, version folded
into the seed. Statistical adequacy is checked empirically in T5;
counter-mode Philox is admitted only if T5 finds bias.

Key = `(version, seed, stream, invocation, logical index, projection
slot)`. Logical index is the element's index in the user's array origin
(subarray calls adjust origins — same element, same word). Projection
slot distinguishes multiple projections inside one logical output
(`ConvertToBlockMaxAbsFinite`'s scale vs elements).

**Projection slots are enumerated in a `const SLOTS` table.** `SLOTS`
lives in `catalogue.jl` alongside operation metadata and maps each
operation to the ordered list of `UInt32` slot labels its element-side
kernel emits per logical output. This is the authoritative list; the
release manifest records a SHA-256 over `SLOTS` so an inadvertent slot
renumbering shows up as a signature mismatch at plan preparation and as
a T5 golden-vector failure:

```julia
# catalogue.jl
const SLOT_SCALE           = UInt32(0)
const SLOT_ELEMENT         = UInt32(1)
const SLOT_TIE_BREAK       = UInt32(2)   # StochasticC ties (§5.7.4)
# ... one const per named slot ...
const SLOTS = Base.ImmutableDict(
    :Convert                       => (SLOT_ELEMENT,),
    :ConvertToBlockMaxAbsFinite    => (SLOT_SCALE, SLOT_ELEMENT),
    :Add                           => (SLOT_ELEMENT,),
    :FMA                           => (SLOT_ELEMENT,),
    :BlockReduceAdd                => (SLOT_ELEMENT,),        # one projection at close
    :BlockDotProduct               => (SLOT_ELEMENT,),
    # ... every operation, every slot it uses ...
)
const SLOTS_SHA256 = "…"   # recorded in the release manifest
```

Callers must not fabricate slot values. Passing an unregistered `UInt32`
to `entropy_word` in a debug build triggers an assertion; the release
build trusts the type-parameter path that only ever forms slots from
`SLOTS[op]`. Adding a slot is a catalogue change: new `const`, new
`SLOTS` entry, updated `SLOTS_SHA256`, T5 golden-vector row per new
slot. This closes the "multiple projections inside one output" clause
above: without an enumerated table, the projection-slot half of the
reproducibility key is folklore.

`SequentialEntropy(rng::AbstractRNG)` consumes in sequence and **forces
the sequential generic plan** — the planner refuses to thread or tile it.
One word per logical output still holds.

Requesting a stochastic mode without an entropy source throws
`EntropyRequiredError` at plan preparation, never mid-loop.

### 12.2 The reproducibility gate (Stage 8 exit)

Same signature, entropy identity, inputs, logical origin ⇒ identical
bytes across: scalar loop; 1 and many threads; different chunk sizes;
dense and packed storage; SIMD and generic loops; subarray calls with
adjusted origins. T5's backbone. This property — schedule-independent
threaded stochastic rounding — is the design's reason for indexed
entropy; an RNG passed into kernels cannot have it.

---

## 13. Blocks, reductions, scaled, packed — `blocks.jl`, `packed.jl`

### 13.1 Exact reductions ([§6.3](IEEE_D1.md#L2244))

Terms: decoded scale×element products — exact `Datum`s. Significand
width depends on op: `BlockReduceAdd/Subtract` and the `Scaled*` family
feed single-scaled products, ≤ 16-bit sigs (assert); `BlockDotProduct`
feeds four-factor products, ≤ 32-bit sigs (assert). Product exponent
spread across the lattice ≤ 506 binades.

The draft defines reductions as left folds of exact ω-operations with a
single final projection. Because the shipped accumulator is **exact**, any
association order produces the same mathematical sum, so partition and
merge order cannot change the result code — the "reassociation is
approximate unless proven" concern discharges by exactness, once, here,
rather than per-optimization.

**Shipping design: `LimbAcc` as the default, windowed as opt-in.**

- **Shipped default (Stage 9):** derived fixed-limb signed
  superaccumulator `LimbAcc{L}` in `src/blocks.jl`. Exact by
  construction; threading-safe by exact integer addition; carries every
  reduction on every input distribution. This is the reduction path
  users see out of the box.
- **Reference (Stage 9, T6-only):** unbounded `BigInt`-limb accumulator
  in `reference/block_oracle.jl`, in the same fixed-point frame as
  `LimbAcc`. Never reachable from production; used solely by T6 to
  check the shipping accumulator (§17.4).
- **Optional windowed adapter (Stage 11, gated):** a serial windowed
  accumulator with bounded restart, faster than `LimbAcc` on
  well-behaved operand distributions. Ships only after clearing three
  independent gates — derived restart bound, exhaustive equivalence,
  and a measured end-to-end win — detailed in §13.1.4 below. If any
  gate fails, the adapter is deleted and `LimbAcc` remains the only
  reduction path.

This ordering resolves the tension between speed and provability:
shipping the exact accumulator first means every user gets a
schedule-independent, threading-safe result on day one, and every future
optimization proves itself against a live oracle rather than a paper
design. Threaded reductions never use the windowed adapter — the
exactness argument for schedule-independence lives with `LimbAcc` alone
(invariant 7, §2). Merging windowed partials across chunks would
reintroduce cross-chunk tail reasoning and forfeit that invariant.

**Shipped accumulator — fixed-limb signed superaccumulator.**

```julia
struct LimbAcc{L}
    pos::NTuple{L, UInt64}   # non-negative running sum, fixed-point at 2^emin_summand
    neg::NTuple{L, UInt64}   # non-negative subtractor; final value = pos - neg
end
```

Limb count derived per op at signature preparation from the **summand**
exponent range — the exponent of what the stream actually adds, not of
the raw operand — with per-op derivations:

```text
emax_summand = emax_scale + emax_element               # BlockReduceAdd/Subtract
emin_summand = emin_scale + emin_element
emax_summand = emax_scale_a + emax_element_a
             + emax_scale_b + emax_element_b            # BlockDotProduct
emin_summand = emin_scale_a + emin_element_a
             + emin_scale_b + emin_element_b

sig_bits   = SIG_BITS_SINGLE_SCALED   # 16 for BlockReduceAdd/Sub
             or SIG_BITS_DOT           # 32 for BlockDotProduct (four-factor product)
span_bits  = (emax_summand - emin_summand) + sig_bits + ceil(log2 B) + 2
L          = ceil(span_bits / 64)
```

The two constants are declared at the top of `blocks.jl`:

```julia
const SIG_BITS_SINGLE_SCALED = 16   # decoded datum × decoded scale carrier
const SIG_BITS_DOT           = 32   # four-factor Datum product (§13.1 dot path)
```

The two constants have distinct physical origins — the single-scaled
summand is a 8-bit-sig datum times an 8-bit-sig scale, giving 16; the
dot summand is a four-factor product of 8-bit-sig datums, giving 32 —
and keeping them separate documents that origin in the source. The
plan-preparation code refuses a mismatched `sig_bits` at G3.

For the byte-domain lattice at `B ≤ 4096`, `L ≤ 9` for `BlockReduceAdd`/
`BlockReduceSub` (using `SIG_BITS_SINGLE_SCALED`) and `L ≤ 17` for
`BlockDotProduct` (using `SIG_BITS_DOT`). `L` is re-derived per op;
the doubled scale-exponent span in `BlockDotProduct` is handled by the
four-factor `emax_summand` formula, not by inheriting from
`BlockReduceAdd`. The earlier draft's formula `(2·emax_datum −
2·emin_datum) + 16` conflated summand and format ranges and was wrong
for `BlockDotProduct`; this replaces it.

**Per-element update.** Each incoming summand is an exact `Datum`.
Position its `sig` at bit offset `off = d.exp - emin_summand` in the
`L`-limb fixed-point frame; `off ≥ 0` by construction of `emin_summand`.
`LimbAcc` is immutable, so `update` returns a new value — the tuple
sits in registers at `L ≤ L_static`, so this costs no allocation:

```julia
@inline function update(acc::LimbAcc{L}, d::Datum) where {L}
    off      = d.exp - emin_summand            # ≥ 0 by construction
    hi_limb  = (off >> 6) + 1                   # 1-based limb index of the LSB
    intra    = off & 0x3f                       # bit offset within that limb
    # Widen d.sig into two limbs so any intra > 0 spans a limb boundary.
    lo_word  = UInt64(d.sig) << intra
    hi_word  = intra == 0 ? UInt64(0) : UInt64(d.sig) >> (64 - intra)
    addend   = _sparse_addend(Val(L), hi_limb, lo_word, hi_word)  # NTuple{L,UInt64}
    if d.neg
        return LimbAcc{L}(acc.pos, add_limbs(acc.neg, addend))
    else
        return LimbAcc{L}(add_limbs(acc.pos, addend), acc.neg)
    end
end
```

`_sparse_addend` is a `@generated` helper producing the `L`-tuple with
zeros in every limb except `hi_limb` (which receives `lo_word`) and
`hi_limb + 1` (which receives `hi_word`, guarded by the `L` bound).
Both halves live inside the summand range by construction of `L` above,
so `_sparse_addend` asserts `hi_limb + 1 ≤ L` in debug builds.

**Register residency and static ceiling.** `NTuple{L, UInt64}` values
decompose into scalar SSA registers under Julia 1.12's SROA reliably only
for small `L` — [large tuples spill and allocate](https://discourse.julialang.org/t/why-do-large-ntuples-allocate/73145).
`L_static = 4` (a `const` in `src/blocks.jl`, overridable via
`Preferences.jl` at build time by callers who benchmark higher and pin
`@allocated == 0` themselves). At `L ≤ L_static` the accumulator, its
per-element updates, and the merge are all on-stack and register-
resident; at `L > L_static` the accumulator lives in a caller-owned
`Vector{UInt64}` workspace of length `2L` (pos and neg concatenated),
reported by `plan_traits` before execution (`WorkspaceError` if
undersized). No silent per-element allocation exists at any `L`.

**Multi-limb addition — carry-propagating, not elementwise.** Tuples do
not support `+`, and broadcast `.+` on `NTuple{L, UInt64}` drops
carries. A carry-propagating primitive is required, but the shape of
that primitive matters as much as its arithmetic.

The naive `ntuple(Val(L)) do i ...; carry = ...; end` closure form is
**wrong**: closing over a mutated `carry` boxes it into a `Ref{UInt64}`
under Julia 1.12's captured-variable rules ([performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured)),
violating W3 and blocking the `adc`-chain lowering that motivates the
primitive. Two idiomatic 1.12 forms produce a straight-line unrolled
`adc` chain instead:

```julia
# Recursive tail-call form; each call is @inline'd and the recursion is
# unrolled at compile time because L is a type parameter.
@inline add_limbs(a::NTuple{L, UInt64}, b::NTuple{L, UInt64}) where {L} =
    _add_limbs(a, b, UInt64(0), Val(1), Val(L))

@inline function _add_limbs(a, b, carry::UInt64, ::Val{I}, ::Val{L}) where {I, L}
    s,  c1 = Base.add_with_overflow(a[I], b[I])
    s2, c2 = Base.add_with_overflow(s, carry)
    ncarry = UInt64(c1) + UInt64(c2)
    if I == L
        return (s2,)
    else
        return (s2, _add_limbs(a, b, ncarry, Val(I + 1), Val(L))...)
    end
end
```

Or, equivalently, a `@generated` function that emits the `L`-limb
straight-line body directly:

```julia
@generated function add_limbs(a::NTuple{L, UInt64}, b::NTuple{L, UInt64}) where {L}
    exprs = Expr[]
    limbs = Symbol[]
    push!(exprs, :(carry = UInt64(0)))
    for i in 1:L
        si  = Symbol("s",  i)
        s2i = Symbol("s2", i)
        push!(exprs, :(($si,  c1) = Base.add_with_overflow(a[$i], b[$i])))
        push!(exprs, :(($s2i, c2) = Base.add_with_overflow($si, carry)))
        push!(exprs, :(carry = UInt64(c1) + UInt64(c2)))
        push!(limbs, s2i)
    end
    push!(exprs, Expr(:tuple, limbs...))
    Expr(:block, exprs...)
end
```

Both forms lower to `L` `adc` instructions on x86-64 under Julia 1.12's
LLVM backend, allocate zero, and share code between per-element updates
and cross-thread merges. Pick one during Stage 9 based on measured
compile time — the recursive form composes with `@code_llvm` more
cleanly, the `@generated` form is easier to review. Do **not** ship the
mutating-closure form.

**Threading is bit-deterministic by construction.** Per-thread
`LimbAcc{L}` accumulators merge by two exact `add_limbs` calls (one for
`pos`, one for `neg`), then a canonicalization step:

```julia
@inline function merge(a::LimbAcc{L}, b::LimbAcc{L}) where {L}
    pos = add_limbs(a.pos, b.pos)
    neg = add_limbs(a.neg, b.neg)
    # Canonicalize: exactly one of {pos, neg} is nonzero after merge.
    # cmp_limbs is an integer-limb comparison, associative and
    # commutative in its output ordering; the branch below is a total
    # function of (pos, neg) and preserves the invariant "|pos - neg|
    # in one lane, zero in the other, sign carried by which lane is
    # nonzero."
    return cmp_limbs(pos, neg) >= 0 ?
        LimbAcc{L}(sub_limbs(pos, neg), ntuple_zero()) :
        LimbAcc{L}(ntuple_zero(), sub_limbs(neg, pos))
end
```

Without canonicalization the merged accumulator can hold two nonzero
lanes and the derived `L` bound (which counts *summand* range, not
*sum-of-partials* range) would no longer suffice: two partials each near
the worst-case span could push `pos` and `neg` past the bound while the
canonical `|pos − neg|` still fits. Canonicalization is the property
that makes the L derivation in §13.1 valid for both the serial and the
threaded path.

Integer addition and subtraction are associative and commutative under
this canonicalization (an easy induction: the pair-merge is
commutative, and any two associations produce the same `|pos − neg|`
lane), so the merged result is bit-identical to the serial one across
any partition, chunk size, or thread count — no reassociation qualifier
is needed, no scheduler flag changes bits. A `SequentialEntropy`
reduction is not threaded; the planner refuses the partition (§12.1).
Should the Stage-11 windowed adapter ship, it remains **serial-only**;
merging windowed partials across chunks is forbidden.

**Reference — `BigInt` oracle in `reference/`.**
`reference/block_oracle.jl` uses an unbounded `BigInt`-limb accumulator
with the same fixed-point frame as `LimbAcc`. T6 checks the shipping
accumulator against it on the full lattice, on cancellation towers
spanning the exponent range, and at every limb-width transition
`L ∈ 1:L_max_dot`, where the two block-operation families have distinct
limb-width ceilings derived from their summand ranges:

- `L_max_sum = 9` for `BlockReduceAdd`/`BlockReduceSub` (summands are
  `Datum` with ≤ 16-bit significand; §13.1 formulas give a span of
  `≤ 549` bits, so `⌈549/64⌉ = 9`).
- `L_max_dot = 17` for `BlockDotProduct` (four-factor `Datum` products
  reach a 32-bit significand and the summand span reaches 1 058 bits,
  so `⌈ 1 058 / 64⌉ = 17`).

T6 therefore covers `L ∈ 1:17` (the union), and every row records which
sub-range applies. Should the Stage-11 windowed adapter ship, T6
additionally checks it against `LimbAcc` on the same input classes.

Streaming shape in both forms: decode scale once; per element decode,
scale, update the accumulator; project once. No decoded temporary vector
is ever allocated.

`BlockReduceAdd` seeds 0; `BlockDotProduct` forms exact per-lane
four-factor products (`Datum` with ≤ 32-bit significand — one limb), then
feeds `LimbAcc{L}` with `L` re-derived per the four-factor `emax_summand`
formula.

**`BlockReduceMultiply` — not a sum.** A product of `B` 8-bit
significands needs up to `8·B` bits and the exponent is a sum of `B`
exponents. It is not a two-term window problem and does not use
`LimbAcc`. Ship it as a separate exact accumulator:

```julia
struct MulAcc{Lm}
    sig::NTuple{Lm, UInt64}   # unsigned exact significand product
    exp::Int64                # exact exponent sum
    neg::Bool                 # parity of operand signs
    saw_zero::Bool            # true collapses to exact zero at final projection
    saw_nan::Bool             # NaN operand seen; dominates on merge
    saw_inf::Bool             # signed-Inf operand seen; combined with signs and saw_zero at projection
end

# Per-element step: multiply the current Lm-limb significand by an 8-bit
# operand, propagate carries, and accumulate exponents. MulAcc is immutable;
# each step returns a new MulAcc. The chain is straight-line @generated on Lm.
@inline function mul_step(acc::MulAcc{Lm}, d::Datum) where {Lm}
    d.sig == 0        && return MulAcc{Lm}(acc.sig, acc.exp, acc.neg, true, acc.saw_nan, acc.saw_inf)
    isnan_datum(d)    && return MulAcc{Lm}(acc.sig, acc.exp, acc.neg, acc.saw_zero, true, acc.saw_inf)
    isinf_datum(d)    && return MulAcc{Lm}(acc.sig, acc.exp, acc.neg, acc.saw_zero, acc.saw_nan, true)
    new_sig = _widen_mul_u8(acc.sig, UInt8(d.sig))       # Lm-limb × UInt8 → Lm-limb
    return MulAcc{Lm}(new_sig, acc.exp + d.exp, xor(acc.neg, d.neg),
                      acc.saw_zero, acc.saw_nan, acc.saw_inf)
end

# Threaded merge: multiply two Lm-limb significands; combine flags.
@inline function merge(a::MulAcc{Lm}, b::MulAcc{Lm}) where {Lm}
    MulAcc{Lm}(_mul_limbs(a.sig, b.sig), a.exp + b.exp, xor(a.neg, b.neg),
               a.saw_zero | b.saw_zero, a.saw_nan | b.saw_nan, a.saw_inf | b.saw_inf)
end
```

`_widen_mul_u8` (per-element, `Lm`-word by `UInt8`) and `_mul_limbs`
(merge, `Lm`-word by `Lm`-word) are both `@generated` on `Lm`, emitting
the `mul`/`adc` chain directly — the exact-multiplication analogues of
`add_limbs`. Both are Stage-9 deliverables gated by W2 on
`BlockReduceMultiply` inputs. T6 rows for multiply cover `saw_zero`,
`saw_nan`, `saw_inf`, signed-parity, cross-limb carry, and thread-merge
commutativity.

`Lm = ceil((8·B_max) / 64)`: `Lm = 1` for `B ≤ 8`, `Lm = 2` for
`B ≤ 16`, etc. Products grow monotonically; each step is an `Lm`-limb
multiply-by-8-bits, cheap and exact. Zero collapses to `saw_zero = true`
and short-circuits the final significand. NaN operands set
`saw_nan = true` and short-circuit the projection to the format NaN.
Signed infinities set `saw_inf = true`; the projection combines
`saw_nan`, `saw_inf`, `saw_zero`, and `neg` under the IEEE-754
special-value table (NaN dominates Inf × Zero). Threaded merge combines
two `MulAcc{Lm}` values by (i) OR-merging `saw_nan`, `saw_zero`,
`saw_inf`; (ii) XOR-merging `neg`; (iii) `Lm`-limb multiply of the two
partial significands with carry propagation; (iv) exponent add. Each
primitive is commutative and associative, so the merge is byte-identical
across partitions and thread counts. `Lm ≤ L_static` shares the static
path; `Lm > L_static` uses the same workspace mechanism as `LimbAcc`.

### 13.1.4 Optional windowed accumulator — Stage 11, gated

A windowed reduction with bounded restart can be faster than `LimbAcc`
serially when the operand distribution is well-behaved (few restarts,
small `L`). It is a Stage-11 opt-in adapter, never a Stage-9
deliverable, and it ships only after clearing three independent gates.
Failure of any gate deletes the adapter; `LimbAcc` remains the only
reduction path.

**(a) Derived restart bound (never a flat count).** Each restart lowers
the anchor by ≥ `63 − hb − 17` binades, giving a bound of
`⌈506/(63 − hb − 17)⌉ + 1`, tabulated across the covered block sizes:

| `B` range | `hb = ⌈log₂ B⌉ + 1` | `63 − hb − 17` | restart bound |
|:---|:---|:---|:---|
| `B ≤ 2` | 2 | 44 | 13 |
| `B ≤ 8` | 4 | 42 | 14 |
| `B ≤ 32` | 6 | 40 | 14 |
| `B ≤ 64` | 7 | 39 | 14 |
| `B ≤ 256` | 9 | 37 | 15 |
| `B ≤ 4096` | 13 | 33 | 17 |

Small **fixed** term counts admit the far tighter §4.4 counting bound
(2 for FAA). A flat bound like "at most 32 restarts" would satisfy the
test and hide the actual property; the derivation is what the gate
checks.

**(b) Exhaustive equivalence with `LimbAcc`.** The Stage-11 acceptance
suite runs the windowed adapter against `LimbAcc` on:

- adversarial cancellation towers spanning the exponent range;
- the full operand lattice for `B ∈ {1, 2, 8, 32, 64, 4096}`;
- every representable `Datum` at each limb-width transition
  `L ∈ 1:L_max`;
- specials, saturation edges, and the parity/entropy grid covered by
  T2/T5.

Every input class produces bit-identical codes. Any discrepancy — even
one — deletes the adapter. This is stronger than T6's tier-agreement
check because it runs before the adapter reaches production, not after.

**(c) Measured end-to-end win.** The Stage-11 benchmark harness
(§18) reports geometric-mean speedup over `LimbAcc` on representative
workloads and worst-case slowdown on the acceptance suite. Both must
clear the harness threshold recorded in `benchmark/thresholds.md`. "It's
theoretically faster" without a certified measurement does not admit
the adapter.

**Windowed pseudocode (Stage-11 reference).** The algorithm below is
frozen here so that the acceptance suite has a target to check against.
`LimbAcc` semantics are unaffected — this appears in `src/blocks.jl`
only if all three gates above pass.

```text
pass 1: max topexp M over nonzero terms
anchor: place bit M at window bit 63 − hb, hb = ceil(log2(B)) + 1
        (headroom so B same-sign additions cannot overflow either
        sign-lane; wpos and wneg are separate)
pass 2: accumulate wpos, wneg :: UInt64 (separate sign sums);
        terms shifted below bit 0 are dropped; track max dropped topexp
        (t_drop) and count
merge:  # UInt64 subtraction wraps; use sub_with_overflow (or the signed
        # form; both are one machine instruction) to carry sign explicitly.
        # diff, borrow = wpos - wneg (native); when borrow, magnitude =
        # -diff (mod 2^64) = ~diff + 1.
        # The debug assert below reads as "the wraparound case where the
        # 64-bit difference is exactly zero AND we borrowed is
        # unreachable": wpos < wneg with diff == 0 mod 2^64 would require
        # wneg = wpos + 2^64, but wpos, wneg are UInt64 values (< 2^64),
        # so the sum wpos + 2^64 is not representable and no such wneg
        # exists. The old inline comment attributed unreachability to a
        # "[0, 2^63] bound" — that was misleading; UInt64
        # REPRESENTABILITY is what makes it unreachable, not any tighter
        # numeric bound the adapter maintains.
        diff, borrow = Base.sub_with_overflow(Int64(wpos), Int64(wneg))
                                              # signed form is equally cheap
                                              # and reads unambiguously
        @assert !(borrow && diff == 0)   # unreachable per the paragraph above
        result_neg   = borrow
        w            = borrow ? (~diff + one(UInt64)) : diff
if no drops                -> exact: project (w, result_neg, anchor exp)
elseif w ≥ 2^(hb + 17)     -> tails cannot reach the rounding decision:
        combined tail sign resolved as in §4.4 — one extra pass over the
        dropped terms in a SECOND window anchored at t_drop (itself
        exact-or-tailed; its drops recurse)
else                       -> RESTART anchored at max(topexp(w'), t_drop),
        w' the exact retained residual
```

Why `sub_with_overflow`: `wpos, wneg ∈ [0, 2^63]` in the adapter's
accounting, so their true difference is in `[-2^63, 2^63]`, which fits
neither `UInt64` (negatives wrap) nor `Int64` (the positive endpoint at
`2^63` overflows exactly). `sub_with_overflow` does the native
subtraction and returns the borrow (or overflow) bit; the sign is
carried in `result_neg` and the magnitude in `UInt64`. Reading the
signed form is enough at the merge site — the borrow flag from the
signed instruction and the unsigned instruction agree in this range,
and the signed form is closer to what a human reads.

Gate **W2** in §17.2 asserts the derived restart bound for the `B` under
test and agreement with `LimbAcc` on constructed cancellation towers,
for every Stage-11 build that includes this adapter. Builds without it
do not run W2's windowed-specific rows.

### 13.2 Block storage

```julia
struct BlockBytes{SF,EF}
    scales::Vector{SF}
    elements::Vector{EF}      # contiguous, block-major lanes
    blocksize::Int
end
```

Constructor checks length divisibility (`DimensionMismatch`).
Structure-of-arrays over bytes; `NTuple{B,EF}` only for small static `B`
and only behind a measured win — SmallFloats measured the tuple-length
allocation cliff (48 B at `B = 8` → 6128 B at `B = 64`); do not
rediscover it. No static-array dependency in the core interface.

### 13.3 Scaled ops and block conversions

`ScaledOp` = `BlockOp` at `B = 1` ([§6.5](IEEE_D1.md#L2357)) — one
definition, so two versions of scaled arithmetic cannot drift; an
optimized scalar adapter is admitted only after exhaustive equivalence.
Scaling is integer significand multiply + exponent add — no wide carrier;
form the complete expression, project once. `ω-BlockDecode` **multiplies**
by the scale; `ω-BlockProject` **divides** ([§6.1](IEEE_D1.md#L2089)) via
`window_div` before projection; scale edge rows (`S = 0 → 0`,
`S = ±Inf → sgn(X)·sgn(S)`, NaN rows) are code-point rows first.
`ConvertToBlockMaxAbsFinite` ([§6.2.3](IEEE_D1.md#L2200)) carries two ρs
(`ρ_s` scale, `ρ` elements) and two projection slots (§12.1). `BlockOp`'s
result scale is an **input**, never invented ([§6.4](IEEE_D1.md#L2310)).

### 13.4 Packed storage

```julia
struct PackedBytes{K}
    words::Vector{UInt64}
    len::Int
end
```

`K = 8` is an identity byte view. Unused tail bits zeroed
deterministically so packed equality and hashing are stable. Fuse
unpack → execute → repack in one traversal. Honest economics: pack
0.53 ns/el, unpack 0.44 ns/el vs a 0.13 ns/el gather — packing pays only
memory-bound. Direct packed-word arithmetic: exhaustive equivalence plus
a measured **end-to-end** win, or it doesn't ship.

---

## 14. Julia surface and external formats

### 14.1 `julia_adapters.jl` — the deliberate `AbstractFloat` surface (Stage 5)

`zero`, `one`, `typemin`, `typemax`, `floatmin`, `floatmax`, `eps`,
`ldexp`; `isnan`, `isinf`, `isfinite`, `issubnormal`, `signbit`, `iszero`;
`==`, `isless`, `hash` (on `(format, code)`, NaN canonical); conversions
to/from Float16/32/64, `BigFloat`, integers; `parse`/`tryparse`; `rand`
(uniform over finite datums — documented choice); `broadcast`/`similar`
with concrete eltypes; zero-copy `reinterpret(UInt8, A)`.

`isless` follows the P3109 total order — **NaN sorts first**, inverting
the Float64 convention; documented loudly, and deliberately: one order for
predicates, sorting, and selection, from one key table. `typemin` on
extended signed is `−Inf`, on finite `M_lo`; NaN is never
`typemin`/`typemax`. `eps(T)` follows the Float64 convention (at one);
`eps(x)` from the ulp at `x`.

### 14.2 Promotion and operators

Promotion targets **Float64** — exact for every datum; a cross-format
`FloatByte` join is not pretended canonical. Base operators require
same-format operands under the immutable package default (NTE,
`SatNone`); cross-format `+` is a `MethodError` with a hint pointing at
`apply`/`Convert`. Scoped override:

```julia
const PROJECTION_SCOPE = ScopedValue{Projection}(DEFAULT_PROJECTION)
with_projection(f, ρ) = with(f, PROJECTION_SCOPE => ρ)
```

Documented as outside the performance interface — explicit
`apply(sig, x, y)` is the fast path; only the Base-operator adapters read
the scope. Benchmark-gate the read; the fallback is deleting the runtime
override. No mutable global numerical policy exists anywhere.

### 14.3 External formats — `external.jl`

An internal `ExternalFormat` adapter trait: decode, encode, format
queries, and a **declaration of which projection modes it implements**. A
host conversion under ambient defaults is not an implementation of
arbitrary `ρ`: the initial conformance path is Float32 under
`(NearestTiesToEven, SatNone)`, implemented and tested at bit level;
other modes are bit-level implementations or explicit unsupported errors.

Bridge semantics ([§5.8](IEEE_D1.md#L865)): external decode collapses
every host NaN to the single internal NaN and both host zero signs to the
single zero; external encode emits positive zero and a quiet NaN with
zero payload. `reinterpret` drives field-level tests; the semantic oracle
is exact integer/dyadic decoding — a `BigFloat` cast may carry finite
values but must never decide NaN/infinity/negative-zero policy.

Float64 and Float32 ship first (the decode atlas already proves datum
exactness in both). BFloat16 is `ext/FloatBytesBFloat16Ext` — no core
dependency for an optional interop type.

---

## 15. `reference/` — the independent authority

Four files, **zero imports from `src/`** (grep-enforced):

- `draft_model.jl` — formats, decode, encode, specials, ordering, written
  directly from draft clauses with line citations; `Rational{BigInt}`
  values plus symbolic specials.
- `rational_projection.jl` — projection by **neighbor search**: enumerate
  the target's sorted finite datums (cached per format), bracket the
  exact rational, decide by §7.2's exact quotient/remainder rules
  (stochastic thresholds as exact `Rational` comparisons against
  `R/2^N`), then the draft's saturation table.
- `interval_oracle.jl` — directed MPFR enclosures with sidedness flags
  for the elementary functions; used by T3 and table certification.
- `block_oracle.jl` — the `BigInt`-limb reduction oracle (§13.1), in the
  same fixed-point frame as the shipping `LimbAcc{L}` accumulator; used
  only by T6. Never reachable from production.

Oracle discipline: reference and production may share trait formulas and
special-code constants; they must not share rounding-decision code or
window code (§7.2). Fixture tables from the draft's annexes are parsed as
reviewed test data — expected results are never generated by the code
under test.

SmallFloats (`K ∈ 3:8`) is the second differential adapter, test-only,
version-pinned. Disagreements resolve against the draft text — never by
majority; each resolution is recorded in
`test/assurance/resolutions.md`. No code is copied from SmallFloats into
`src/` or `reference/`.

---

## 16. Conformance registry — `conformance.jl`

Conformance is about operation **specializations**, not exported names.
Immutable data:

```julia
struct Specialization
    operation::DataType
    operand_formats::Tuple
    result_format::DataType
    projection::Projection
end

abstract type Accuracy end
struct ExactAccuracy  <: Accuracy end
struct KappaSteps     <: Accuracy; steps::Int; end
struct KappaInfinity  <: Accuracy end
struct KappaNaN       <: Accuracy end
```

Narrow public query surface — three named predicates, each with a
committed meaning:

- `conforms_exact(spec)::Bool` — an **exact** implementation of `spec`
  is registered and its verification gate has passed. Returns `false`
  for any approximate implementation, even one with `κ = 0` measured
  on the covered set. This is the predicate library code uses when it
  needs behavioral guarantees.
- `supports_any(spec)::Bool` — some implementation of `spec` is
  registered (exact **or** approximate). Use this for discovery,
  documentation, and "is this operation available at all" checks.
- `implementation(spec)` — returns `nothing` or an immutable
  declaration containing implementation identifier, accuracy
  (`ExactAccuracy`, `KappaSteps`, `KappaInfinity`, `KappaNaN`),
  package version, pinned draft revision, and optional attestation
  artifact. Callers dispatch on the declaration's accuracy field for
  finer-grained decisions than the two Boolean predicates.
- `conformance_report()` — the serializable, stable-order dump used in
  release manifests and Stage-9 gates.

The old single `supports(spec)` name is retired. Its meaning was
ambiguous — exact-only or any-registered? — and downstream code (T4,
Stage-9 gate, `approximation.jl`) was inconsistent about which it
wanted. `conforms_exact` and `supports_any` name both variants
explicitly; a Stage-9 grep gate refuses any occurrence of the bare
`supports(` identifier in `src/`.

Rules:

- **Never infer exactness from the presence of a method.** Registration
  happens only after the specialization's verification gate passes.
- Approximate implementations carry identifiers distinct from exact
  operations (`approximation.jl`), and no approximate implementation is
  reachable from a defined plan — a type-level test, not a convention.
- Finite κ uses result-format **order-key** distance, never raw code
  subtraction across the signed split. Ladder: NaN agreement, then
  infinity agreement (`MatchOnInfinity`), then max finite value-step
  distance; multi-result rule NaN ≻ ∞ ≻ max ([§5.4](IEEE_D1.md#L468)).
- Every κ or conformance claim states its scope: exhaustive, proved, or
  sampled — a certificate type distinguishes them.
- Human-readable conformance docs are generated from the registry, so
  code, tests, and claims share one source of truth.

Required conforming core ([§5.5](IEEE_D1.md#L521)), expanded by the
registry generator into concrete tuples under `(NearestTiesToEven,
SatNone)`:

```text
F4 = {Binary4p2sf}    F8 = {Binary8p4se, Binary8p3se}
FX = {Float32} initially    Fs = {Binary8p1uf} for scaled ops
```

Families: conversion across the union; negate/abs; reciprocal; add,
subtract, multiply; FMA and FAA; all ten min/max variants; five
comparisons; the predicates and both neighbour ops; all format-level and
extremal queries; scaled add/subtract/multiply with `Fs` scale. Before
release, diff the generated list against the pinned draft's normative
clause — the prose checklist aids implementation and is not the
declaration.

---

## 17. Verification

### 17.1 Tier map → concrete test entry points

| tier | entry | contents |
|:---|:---|:---|
| T1 | `test/lattice/` | all 13 296 codes: `fromcode`/`codepoint` round-trip, encode∘decode bijection both directions, decode == reference (Float32 and Float64 exactness), class partition (every code exactly one class), order keys monotone vs reference, neighbour graph laws, parse/show round-trip, signed symmetry around the NaN code |
| T2 | `test/projection/` | per format × ρ: every grid edge, midpoint, midpoint±1 unit, saturation edge, exact point, special kind; the full §7.3 saturation cross-product; each rounding mode; the `P == 1` parity helper; selected entropy words incl. `0` and `2^N − 1` |
| T3 | `test/semantics/` | exhaustive unary (all codes × ops × formats); exhaustive same-format binary (**2 504 832 pairs**, footnote below) vs rational reference; FAA associativity; cancellation and exponent-gap suites; RSqrt/Hypot exhaustive per format; interval-checked elementary |
| T4 | `test/plans/` | every materialized table entry == `defined_code`; every plan class byte-identical on identical inputs; planner cost-only assertion |
| T5 | `test/stochastic/` | exhaust all `R` at small `N` (fractions below/at/above half); max-`N` edge pins incl. the `N > 32` path; §12.2 reproducibility matrix; entropy-word golden vectors |
| T6 | `test/blocks/` | shipping `LimbAcc{L}` vs `reference/block_oracle.jl` (`BigInt`) on cancellation towers, exponent extremes, specials, `B ∈ {1, 2, 8, 32, 64, 4096}`, thread partitions, and every limb-width transition `L ∈ 1:17` (union of `L_max_sum = 9` for BlockReduceAdd/Sub and `L_max_dot = 17` for BlockDotProduct; §13.1); scale edge cases; left-fold law. When the Stage-11 windowed adapter ships, T6 adds windowed-vs-`LimbAcc` rows on the same inputs |
| T7 | `test/assurance/failure_injection.jl` | corrupt table digests, undersized workspaces, invalid codes via raw `reinterpret`, forced `InconclusiveProjectionError`, aliasing violations |
| T8 | `test/julia_interface/` | `@inferred` across the surface, allocation counts, promotion, generic algorithms (`sort`, `sum`, `Dict`), hashing, doctests |
| release | `test/assurance/release.jl` | + full 177 M cross-format binary lattice; full same-format ternary (564 M); ~45 min, env-gated (`FLOATBYTES_RELEASE_TEST=1`), sharded by high operand byte in scheduled CI |

**T3 pair-count derivation.** Same-format binary coverage across the
byte domain is `Σ_F |F|²` over the 120 covered formats. The sum
decomposes as: 8-bit formats contribute 60 formats × 256² = 3 932 160
raw entries, 7-bit formats contribute 30 × 128² = 491 520, and formats
≤ 6 bits contribute the balance. After deduplication (codes
canonicalised via `canonical_code` — the invariant-11 view that
collapses signed-`−0` / unsigned-`0` and normalises trap bits) and
after excluding cross-signature pairs that appear in the ternary suite
instead, the effective enumerated set is 2 504 832 pairs. Row-by-row
justification lives in `test/semantics/pair_count.md` alongside the
enumerator. T3's pair-count assertion runs at test start-up and blocks
the suite on any mismatch — the number is a checked invariant, not a
legend.

A requested tier runs in full or errors — never silently downgrades.
Comparison is **code identity in identical logical order**; numerical
tolerance is not evidence for a finite code function. Failure output
includes format, ρ, raw operand codes, decoded operands, reference and
candidate results — a bare decimal mismatch cannot reproduce a rounding
failure.

### 17.2 Gates

- **W1** — window-lemma enumeration as a test: Add pairs over all 120
  formats; FMA triples per ComposedPlan §1.1 **plus the residual gap
  `K ∈ {7,8}, P ≥ 2` plus FAA triples at the same coverage**. Release
  tier; a fixed sampled slice in CI. The plan's evidence scripts are not
  in this repository — Stage 0 recreates `w1_window_lemma_add.jl` and
  `w1_window_lemma_fma.jl` under `test/assurance/` and adds a new
  `w1_window_lemma_faa.jl` with the same `Rational{BigInt}` differential
  structure (all `K ≤ 6` exhaustive, plus the widest-spread `K ∈ {7,8}`
  formats). FAA was not part of the original evidence; it is enumerated
  fresh, not extended by analogy.
- **W2** — Stage 9 requires the shipping `LimbAcc{L}` accumulator to
  agree with the `BigInt` oracle (`reference/block_oracle.jl`) on the
  full T6 input classes; threaded merge is partition-independent across
  every tested thread count and chunk size. If the Stage-11 windowed
  adapter is built, W2 additionally asserts its derived restart bound
  (§13.1.4) and bit-identical agreement with `LimbAcc` on adversarial
  cancellation towers.
- **W3** — zero allocation on every defined path, unconditionally:
  `@allocated == 0` (after warm-up) for every catalogue op at
  representative and adversarial operands, every format class, plus the
  dense loop per element. No "except when" clause exists to test.
  W3 has a **static rider** on top of the runtime check:
  `Base.Compiler.EscapeAnalysis` is consulted on the typed IR of every
  hot-path method, and the rider asserts that the `NTuple{L, UInt64}`
  fields of `LimbAcc{L}` (§13.1) and `MulAcc{Lm}` do not escape their
  update method. The runtime check is a sample; the static check is an
  invariant, and disagreement between them is a Stage-9 gate failure
  — either the runtime harness missed a case or an EA regression in
  the compiler removed an optimization the plan depends on.
- **G9** — traits and queries constant-fold at every alias.
- **G10** — surface totality, wide and shallow: every catalogue op ×
  every legal arity/format combination returns and holds its declared
  result type. Generated from the catalogue so it cannot drift. Kept
  deliberately shallow — its value is totality.

### 17.3 Stage-4 entry condition (risk concentration)

No kernel is trusted before W1 is complete over the residual gap and FAA.
If the window lemma fails on any unreached case, the thesis is revisited
before more code is written — Stage 4 is deliberately small and its gate
deliberately the most expensive relative to the code it covers.

### 17.4 Adapter equivalence chain

```text
reference == defined_code == direct byte law == local table
          == dense == generic == SIMD == threaded == packed
```

One harness walks the chain pairwise on identical inputs in identical
logical order; adding an adapter adds a link, not a harness.

### 17.5 Static checks

Concrete inferred result types on everything exported; zero warm
allocations where promised; no dynamic dispatch, `Any`, lock, dictionary,
or symbol lookup inside loops (code-inspection scans of hot functions);
no `BigInt`/`BigFloat` on any fast path; stochastic state passed as data;
tables bounds-safe and type-stable; bounded table size asserted before
allocation; workspace validated before destination mutation.

### 17.6 Grep-enforced structural rules (one testset)

- `reference/` never imports `src/`.
- `reinterpret(` into `FloatByte` only in `domain.jl` and
  `execution/dense.jl`.
- `BigFloat`/MPFR under `src/` only in `operations/elementary.jl`.
- The string `SmallFloats` never appears under `src/`.
- The load-order comment in `FloatBytes.jl` matches the `include` order.
- Every native-code producer is `saturate∘encode`; no other function in
  `src/` returns a `UInt8` typed as a native code. Raw `reinterpret`
  remains the explicit nonnumeric exception, confined to `domain.jl` and
  `execution/dense.jl` above.
- Approximate-operation type parameters (§16) cannot populate exact
  plan slots — a type-level test, not a convention.
- Generated aliases and artifacts match their generators and recorded
  SHA-256 digests.
- Docstring bodies under `src/` contain no `\(`, `\[`, or `\text`
  escape sequences (§5.6). These render a package unloadable when the
  parser rejects the escape; the grep is scoped to docstring content
  only, not to code strings that legitimately contain backslashes.

### 17.7 Four-axis review

Every stage passes review along four independent axes. A gate green on
correctness but silent on integrity is not a release-quality gate.

**Correctness.** Ordered special cases; exactly one projection per
result; all four format quadrants exercised (signed/unsigned ×
finite/extended); `P == 1` handled by its named helpers; subnormal,
zero, and canonical-byte invariants asserted; every check backed by
independent evidence in `reference/`.

**Performance.** Format decisions are constant-folded (G9); hot paths
are `@inferred` and allocation-free (W3); array/block traversals are
streamed with byte-view access; tables are cache-aware and their build
cost recorded; end-to-end wins are measured on representative workloads,
not only inner instruction counts.

**Integrity.** Canonical construction only (`fromcode` and its
plan-owned neighbors); entropy is explicit and typed; exact and
approximate names are separated at the type level; artifacts are
reproducible and digest-checked; performance contracts are truthful
(median-vs-min separation, allocation counts, cache expectation); no
partial mutation on error paths.

**Clarity.** One owner per behavior; explicit policy over implicit
default; consistent vocabulary across `src/`, `reference/`, and
`test/`; small, review-scoped interfaces; reproducible failures with
complete provenance; objective stage gates from the §19 table — never
"looks good".

---

## 18. Benchmarks — `benchmark/`

### 18.1 Harness rules (paid-for doctrine, adopted whole)

- One variant per process; A/B via a git worktree pinned at the base
  commit, runs in immediate alternation (machine drift between blocks
  reached ~10%; a multi-variant process once reported a 1.87× win that
  process isolation reversed into a 31% loss).
- Preflight aborts the run if a warm scalar path allocates.
- Accumulate into a live value; end with `Base.donotdelete` — a loop that
  overwrites its accumulator measures dead-code elimination.
- No closures over non-`const` globals in timed code (distorts ratios,
  not just absolutes).
- Warm up until numbers stabilize; seed all inputs.
- Latency and throughput reported separately; compile time and code size
  measured as well as steady state.
- Every row names: source table+row, ns/el and cycles/el, elements/s,
  bytes moved, working set, expected cache level, allocations, threads +
  parallel efficiency, CPU model, state (cold/prepared/warm). Mixing
  states is a reporting error.
- Input suites: uniform random codes, sequential, constant, sparse
  zeros, skewed, alternating extrema+specials, cancellation-heavy,
  exact-vs-inexact mixes; contiguous/strided/offset/packed/aliased
  layouts; homogeneous and heterogeneous formats; special-heavy data.

### 18.2 Targets

ComposedPlan §15.1 is the target table, adopted verbatim — each row names
its source. Acceptance rows for this guide's algorithms:

- ternary `vmap` (shipped compute path): ≤ 3 ns/el, **0 allocs**;
- wide-spread 64 K table build: **≤ 80 µs, ≤ 44 allocs** — format-
  independent because no escalation exists;
- scalar medians ≈ minima: the medians are the story — an integer kernel
  has no escalation distribution, so the target is not primarily a lower
  minimum but a median that equals it.

An optimization that misses its target is deleted. Dormant adapters still
cost tests, compilation, and maintainer knowledge. A benchmark win never
becomes a registered exact adapter until exhaustive differential testing
passes.

---

## 19. Stages

Each stage lands green and reviewable as a small PR series; later stages
depend only on earlier interfaces. Optimization begins only after
independent defined behaviour exists. Performance work never lands in the
same commit as a semantic change (the A/B harness depends on the
separation).

| stage | build | exit gate |
|:---|:---|:---|
| 0 | Phase-0/CI policy (§3.3); error taxonomy stubs; `test/stage0_envelope.jl` machine-checking every §1 number; benchmark + verification manifests; recreate + extend W1 evidence scripts | envelope green on all 120 formats; empty suite green on 1.12 + stable; draft revision pinned and visible |
| 1 | `reference/` complete; golden corpus (decode+class+order for all 13 296 codes) with digests; SmallFloats differential adapter | reference agrees with SmallFloats on the corpus; disagreements resolved against draft text and logged |
| 2 | `domain.jl`, `formats.jl`, `show.jl`: type, traits, `fromcode`, atlas, order keys; atlas-vs-tables benchmark decided | T1 exhaustive vs reference; layout + inference pinned; **no arithmetic** |
| 3 | `policies.jl`, `projection.jl`: all modes, all saturation cells, stochastic predicates | T2 exhaustive vs `rational_projection.jl`; Projection is the only encoder |
| 4 | `result_protocol.jl`; seam; spine kernels (Convert, sign ops, Add/Sub, Multiply, Divide, Sqrt/Recip/RSqrt, FMA, FAA, Hypot, Min/Max core, Exp/Log/Sin via §9) | **W1 full (incl. residual gap + FAA)**; **W3**; T3 same-format exhaustive for the spine |
| 5 | `julia_adapters.jl` scalar surface, `external.jl` Float32/Float64 bridges, default ρ, `with_projection` | T8 scalar; bit-level external special-value tests; inference + allocation baselines |
| 6 | unary tables + checksums; three plan classes; dense + generic execution | every entry == reference; build alloc target; T4 |
| 7 | full catalogue: selection, comparisons, predicates, Class, neighbours, format-level ops, remaining elementaries | G10 full width; generated surfaces provably from catalogue |
| 8 | `execution/stochastic.jl`: indexed entropy, decision records, golden vectors | T5 covers §12.2 reproducibility matrix plus stochastic-mode T5 rows (golden entropy vectors, `N > 32` edge case) |
| 9 | `blocks.jl`: `LimbAcc{L}` exact reductions and threaded merge; `MulAcc{Lm}` for `BlockReduceMultiply`; `reference/block_oracle.jl` (BigInt oracle); scaled ops; block conversions | **W2** (`LimbAcc` vs oracle, partition-independence); zero allocation at every tested `B` and every `L, Lm ≤ L_static`; T6; scaled ops verified as `B == 1` blocks |
| 10 | `packed.jl` + fused traversals | byte-identical through the chain; end-to-end target or deleted |
| 11 | one at a time: `BinaryTablePlan`, SIMD, threading, direct-code laws, **optional windowed reduction adapter (§13.1.4)** | each: (a) derivation or specification, (b) exhaustive equivalence with the Stage-9 shipping path, (c) measured end-to-end win — all three or deleted |
| 12 | `approximation.jl` (distinct names), `conformance.jl` registry + κ certificates | no approximate implementation reachable from a defined plan (type-level test); registry reflects tested code only |
| 13 | release review | full matrix (177 M cross-format + full ternary) on x86-64 and AArch64; §16 matrix generated and green; attestation review; deletion test on every seam; a clean checkout reproduces tests, declarations, docs |

Stages 1–7 use SmallFloats as a test-only pinned dependency; once the
exhaustive differentials are green and their digests recorded as goldens,
the dependency is a digest, not a package.

---

## 20. Errors and integrity rules

| error | thrown when | carries |
|:---|:---|:---|
| `InconclusiveProjectionError` | rigorous loop hits its precision ceiling undecided | signature, operand codes, precisions tried, final enclosure; destination unmutated |
| `EntropyRequiredError` | stochastic ρ, no entropy source, at plan preparation | signature, mode |
| `AliasingError` | overlapping views where semantics would differ | array descriptions |
| `WorkspaceError` | caller-owned workspace too small, before any mutation | required vs provided sizes |
| `ArgumentError` | invalid format parameters, out-of-range `fromcode`, invalid stochastic `N` | standard |
| `DimensionMismatch` | block element-count violations | standard |

Mathematical out-of-domain inputs return the specified NaN datum — they
do not throw; there are no flags and no traps. Internal encode
precondition failures are development assertions and indicate package
bugs.

Standing integrity rules (merged from both guides):

- no mutable global numerical policy; no ambient rounding mode or RNG
  dependency;
- no undeclared approximation under an exact operation name;
- no host NaN payload or negative-zero leakage across any bridge;
- no intermediate projection inside fused or reduction operations;
- no conformance entry without tests tied to that exact specialization;
- no generated table without reproducible provenance and reference
  comparison;
- inconclusive is not approximate: no error path writes a guess.

---

## 21. Documentation and release artifacts

User documentation explains: parameter→alias mapping; raw construction vs
mathematical projection; explicit result format and ρ; special-value
differences from IEEE-754 (one NaN, no `−0`, `1/0 = NaN`, NaN-first total
order, `NextGreaterThan(+Inf) = NaN`, saturation as a first-class ρ
component); exact vs approximate identifiers; querying supported
specializations; block scale behavior at zero/infinity/NaN; performance
characteristics and table/code-size tradeoffs.

CI-generated artifacts: the conformance report from the registry;
operation/format coverage tables; benchmark summaries for optimized
adapters; the draft revision + package commit used for attestation;
exhaustive-test manifests and digests for release builds.

---

## 22. Definition of done

FloatBytes 1.0 requires, beyond ComposedPlan §20's fifteen criteria:

- the §16 registry matrix generated with every entry green and attested;
- Float32 interop covered by bit-level special-value tests;
- scaled operations verified as `B == 1` block operations;
- the conformance report naming the pinned draft revision;
- a clean checkout reproducing tests, generated declarations, and docs.

Anything short of these gates can ship as an experimental release, but it
is labeled as such, never described as conforming.

---

## 23. Merge resolutions

Where the two source guides offered alternatives, the selection and its
reason — recorded so neither parent needs consulting to know what
happened.

| topic | window-led guide | phase-led guide | resolution |
|:---|:---|:---|:---|
| value representation | 1-byte `primitive type`, zero-copy byte views | `FloatValue{F,U}` struct, `K ≤ 128` | **primitive type.** ComposedPlan's decided architecture; byte views are the performance backbone; `K > 8` is a non-goal, so the generality buys nothing here. The struct guide's checking spirit lands in `fromcode` being the only safe constructor |
| marker type names | `SignedFmt`/`UnsignedFmt` | `SignedFormat`/`UnsignedFormat`, `FiniteDomain`/`ExtendedDomain` | **phase-led names** — same Base-collision fix, clearer words; invisible to users behind the 120 aliases |
| exact evaluation | 64-bit integer window, proven bounds | `BigInt` Dyadic reference, fixed-width kernels later | **window ships; arbitrary precision is the reference.** Theorem A makes the window sufficient for all fixed 2-term algebra; the phase-led guide's exact quotient/remainder rounding becomes `reference/`'s independently-coded oracle — its best contribution, §7.2/§15 |
| oracle discipline | differential vs `reference/` | "reference and fast must not share rounding code" | **both** — the sharing prohibition is now explicit and grep-adjacent (§7.2, §15) |
| stochastic rounding | mode singletons + external entropy word, indexed | `bits` stored inside the mode struct | **indexed entropy.** Bits-in-mode conflates policy identity with per-call entropy and cannot give schedule-independent threaded results; the phase-led guide's validation rules (`N > 0`, word `< 2^N`, no global RNG) are kept |
| `AbstractFloat` subtyping | subtype now | defer until policy designed | **subtype, staged.** Layout and zero-copy were verified under `AbstractFloat`; the deferral concern is honored by Stage-5 deliberate surface + G10 totality + explicit `MethodError`s before then |
| conformance | thin mention | full registry, accuracy types, attestation | **phase-led registry adopted whole** (§16), harmonized with the catalogue |
| external formats | thin mention | Float32 adapter, NaN/−0 collapse, mode declarations | **phase-led adopted** (§14.3); Float64 added since the atlas already proves exactness |
| file layout | ComposedPlan §16 tree | `internal/` + `fast/` split | **ComposedPlan tree** + `external.jl` and richer `conformance.jl`; reference-vs-fast separation is `reference/` vs `src/`, stronger than an intra-`src/` split |
| stages | 14 ComposedPlan stages | 9 phases with CI-first Phase 0 | **stages, with Phase 0's CI/pinning deliverables folded into Stage 0** and phase exit criteria merged into stage gates |
| Julia floor | 1.11 (`ScopedValue`) | 1.10.10 | **1.12** per project direction; enables `Memory`, `public`, `ScopedValues` without shims |
| reduction order | exact accumulator, order-free | preserve left fold; reassociation ≈ approximate | **exact accumulator** — exactness discharges the reassociation concern once and for all (§13.1); the left-fold law remains a T6 test of the reference |
| reduction implementation | one shipped windowed accumulator + one reference | one exact `BigInt` reference | **`LimbAcc{L}` ships as the default reduction path** in `src/blocks.jl` (exact by construction, threading-safe by integer associativity); the `BigInt` oracle lives in `reference/block_oracle.jl` and is T6-only. A **windowed accumulator is a Stage-11 opt-in adapter** (§13.1.4), admitted only after clearing derivation, exhaustive equivalence, and measured-win gates. Shipping `LimbAcc` first means every user gets schedule-independent results on day one; optimization proves itself against a live oracle |
| reduction → projection handoff | (Plan A) `AlignedX` carrier reused | (Plan B) dedicated `AlignedResidue{Lhi}` carrier + `project_code` method | **`AlignedResidue{Lhi}` adopted** (§4.1). A superaccumulator's below-boundary residue can span a full 64-bit limb, which `AlignedX`'s 16-bit `Datum` remainder would silently truncate. `AlignedResidue` names the carrier at the type level and gets its own verified projector method (§7.1) |
| non-negotiable rules | scattered across sections | (Plan B) numbered list of 12 invariants | **numbered list adopted** at §2; every downstream section cites the invariant it depends on. Concentrating the rules in one place makes PR review mechanical |
| review discipline | tier map + gates | (Plan B) four-axis review (correctness, performance, integrity, clarity) | **both** — the tier map (§17.1) drives what runs; the four-axis review (§17.7) drives what a green tier proves |
| Stage-11 windowed gating | "three gates" phrase | (Plan B) (a) derived restart bound / (b) exhaustive equivalence / (c) measured end-to-end win | **Plan B's (a)(b)(c) framing** at §13.1.4. Each gate is a separately checkable artifact; failure of any one deletes the adapter |
| Julia 1.12+ feature use | mentioned inline | (Plan B) closed-list of language features | **closed-list adopted** at §3.4. Any additional 1.12 feature needs a plan amendment |

Corrections carried from the window-led guide's review (details in §4,
§7): FAA routed through `LimbAcc{3}` (§4.4) with W1 extension; RSqrt
direct integer method with proven at-most-one correction (§4.7); Hypot
enclosure + exact boundary comparison; stochastic `N > 32`
exact-remainder path; derived (not flat) reduction restart bound;
no-`UInt8`-wraparound rule in code-point predicates.

---

## 24. One-page build order

```text
0  prove the numbers          stage0_envelope.jl over all 120 formats; CI + pin
1  build the judge            reference/ (rational + MPFR), goldens, SmallFloats diff
2  build the bytes            FloatByte, traits, fromcode, atlas, order — no math
3  build the ruler            projection: round → saturate → encode, every mode
4  build the window           result_protocol + spine kernels; W1/W3 or stop
5  make it Julia              scalar surface, Float32/64 bridges, default ρ
6  make it fast where free    unary LUTs + three plans + dense byte loops
7  make it complete           full catalogue, G10 wide
8  make chance reproducible   indexed entropy, T5 matrix
9  make many exact            LimbAcc{L} exact reductions + MulAcc{Lm}; BigInt oracle in reference/block_oracle.jl; W2 partition-independence
10 make it dense              packed, fused, or deleted
11 make it faster, one at a time   optional windowed reduction adapter + SIMD + threading + BinaryTablePlan; each: derivation + exhaustive equivalence + measured win, or deleted
12 make claims honest         approximation named apart, registry attested
13 prove it all again         release tier, both architectures, delete the losers
```
