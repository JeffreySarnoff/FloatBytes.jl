# FloatBytes.jl — Implementation Status

Status of the implementation against
[ImplementationPlan.md](ImplementationPlan.md), as of 2026-08-05.
The default test suite runs ~377 000 assertions green on Julia 1.12.6
(`FLOATBYTES_EXHAUSTIVE=1` widens the enumerations).

## Stage coverage

| stage | plan deliverable | status |
|:--|:--|:--|
| 0 | envelope + CI policy + pinned draft revision | **done** — `test/stage0_envelope.jl` machine-checks every §1 number (120 formats, 13 296 code points, bias table, specials, `decode(one) == 1`, `Binary8p1uf` range, Float32/Float64 exactness); `FloatBytes.DRAFT_REVISION` pins the draft copy |
| 1 | `reference/` + goldens + SmallFloats differential | **core done** — `reference/draft_model.jl` (Rational{BigInt} + symbolic decode, exact `√`-value kind) and `reference/rational_projection.jl` (literal §4.7 transcription: RoundToPrecision predicates, saturation rows, encode, ω-op rows), zero imports from `src/`. SmallFloats differential **deferred** (external dependency; the plan makes it test-only) |
| 2 | `domain.jl`, `formats.jl`, `show.jl` | **done** — `FloatByte{K,P,S,D}` primitive type, folding traits, `fromcode`/`codepoint`, format lattice (`format_id` bijection), `Memory`-backed Float64+Float32 decode atlases, generated 120 aliases, non-throwing `show`. T1 exhaustive vs reference. Order-key `UInt16` table deferred: ordering/comparisons run on exact decoded Float64 (provably equivalent; a T1 testset checks the total order) |
| 3 | `policies.jl`, `projection.jl` | **done** — all six deterministic modes, StochasticA/B/C (`N ∈ 1:63` validated), three saturation modes, `round_to_precision → saturate → encode` never merged, `P == 1` code-parity helper (`CodeIsEven`). T2 exhaustive-grid vs the rational reference incl. full small-`N` `R` grids and `N = 32` edges |
| 4 | `result_protocol.jl` + seam + spine kernels | **done** — window substrate (`Datum`, `Aligned`, `window_sum/mul/fma/div/sqrt/rsqrt/hypot`), `LimbAcc{L}` limb accumulator with carry-propagating recursive `add_limbs`, FAA routed through `LimbAcc` (§4.4), all §8.2 kernels with draft-ordered special rows, `defined_code` seam + `apply`. T3 exhaustive small-K binaries, exhaustive K ≤ 5 unaries, sampled 8-bit spot, FMA/FAA triples + associativity |
| 5 | scalar surface + external bridges | **done** — `julia_adapters.jl` (identities, predicates, comparisons, NaN-first `isless`, hash, promotion→Float64, Base operators under `with_projection` scope, parse/rand), `external.jl` (bit-exact Float64/Float32 bridge, −0/NaN collapse, exact Rational/BigFloat/Integer quantize). T8 green incl. `@inferred` and zero-allocation samples |
| 6 | tables + plan classes + dense execution | **not started** (semantics complete without them; tables are cost-only by invariant 7) |
| 7 | full catalogue | **mostly done** — all §8 algebraic, all ten min/max variants, Clamp, comparisons, predicates, Class, neighbours, all 25 unary transcendentals + ArcTan2/ArcTan2Pi. `catalogue.jl` metadata + generated G10 harness deferred |
| 8 | stochastic execution | **core done** — `IndexedEntropy`, splitmix64 `entropy_word` with `Val{N}`, projection-slot constants, golden vectors frozen in T5. Array-level `entropy_for` plumbing awaits Stage-6 dense execution |
| 9–13 | blocks, packed, SIMD/threading, conformance registry, release | **not started** (LimbAcc substrate for Stage 9 already ships in `result_protocol.jl`) |

## Deliberate deviations from the plan text

1. **Stochastic `N > 32` is an explicit error, not implemented.** The
   `AlignedX`/`ExactResidue` exact-residue carriers (§4.1, §4.5, §7.1) are
   deferred; `apply` under `Stochastic*{N > 32}` throws `ArgumentError`
   before any computation ("errors do not write guesses", invariant 12).
   All deterministic modes and stochastic `N ≤ 32` are complete and
   reference-verified.
2. **Sqrt/RSqrt/Hypot windows are wider than the plan's sketches.**
   §4.6/§4.7 as drafted produce ~23–29 significant bits, which
   contradicts §7.1's `k ≥ 33` claim and is insufficient for stochastic
   `N ≤ 32`. The shipped kernels compute at `UInt128` scale
   (`isqrt128`): ≥ 42 significant bits for sqrt/rsqrt, ≥ 49 for hypot.
   Consequence: the open-unit-interval argument makes §4.8's step-5
   boundary refinement unnecessary — no rounding boundary can be interior
   to `(w, w+1)` at that width, so the sticky representative is decisive
   for every implemented mode. Exhaustive T3 agreement with the exact
   rational reference is the checked witness.
3. **`cmp_mag` does not order magnitudes exactly across different
   exponents** (the plan's §4.1 claim is wrong for non-normalized
   `(exp, sig)` pairs). `window_sum`'s exact branch computes the
   symmetric difference instead, so misordering is harmless; the drop
   branch is unreachable under misordering (separation > 24 binades).
4. **FAA uses a fixed `LimbAcc{5}`** (covers the worst K ≤ 8 span,
   `Binary8p1uf`) instead of a per-signature derived `L` — type-stable
   without `Val` plumbing; ≤ 3 spare limbs for small formats.
5. **`LimbAcc` lives in `result_protocol.jl`, not `blocks.jl`** — FAA
   (§4.4) needs it and the substrate must stay dependency-free. Stage 9
   will reuse it from there.
6. **Hypot NaN rows precede the infinity rows**, following the repository
   draft's pattern table literally. (Plan §8.2 claims "Inf rows dominate
   NaN rows here"; the repo draft D1 lists NaN rows first. The draft
   text wins.)
7. **Order-key tables, plan classes, unary result tables** (§5.5, §11)
   deferred with the stages that need them; ordering and selection use
   exact decoded-Float64 comparisons, which T1/T3 prove equivalent.
8. **Structural grep rules** (§17.6) are enforced for the implemented
   subset, with the MPFR rule expressed as "no
   `setprecision`/`setrounding` outside `operations/elementary.jl`"
   (external.jl's BigFloat *conversion surface* is not MPFR machinery).
9. **Transcendental asymptotes**: `Tanh(|x| ≥ 64)` and
   `ExpMinusOne(x ≤ −64)` return the "just inside ±1" representative
   symbolically (the enclosure loop cannot terminate there); all other
   exact points fall out of MPFR's exact-result behavior or explicit
   dyadic rows (`TanPi` poles and quarter points, `ArcSinPi(±1)`, …).
10. **`SLOTS` catalogue table** (§12.1) reduced to the three slot
    constants; the full per-operation table and its SHA-256 arrive with
    the catalogue generator (Stage 7 completion).

## Test suite map (default tier ↔ plan tiers)

| plan tier | file | scope in default run |
|:--|:--|:--|
| Stage-0 | `test/stage0_envelope.jl` | full |
| T1 | `test/lattice.jl` | all 13 296 codes, all laws |
| T2 | `test/projection_t2.jl` | K ≤ 5 full grid + four 8-bit formats; all 6 modes × 3 saturations; stochastic `N ∈ {1,3}` all `R`, three variants, `N = 32` pins |
| T3 | `test/semantics_t3.jl` | K ≤ 4 exhaustive binaries (2 ρ), K ≤ 5 exhaustive unaries (3 ρ), 8-bit sampled, FMA/FAA sampled + special triples + associativity, cross-format Convert, special-row pins, selection family |
| T3 (elementary) | `test/elementary_t3.jl` | three formats × 25 ops × all codes, endpoints re-derived at 512 bits and projected through the reference |
| T5 | `test/stochastic_t5.jl` | entropy golden vectors, contract errors, mean-drift law |
| T8 | `test/julia_interface_t8.jl` | inference, sorting/hash/Dict, promotion, scope, allocations |
| §17.6 | `test/structural.jl` | grep rules, load-order comment |

`FLOATBYTES_EXHAUSTIVE=1` raises T2 to all 120 formats and T3 binaries to
K ≤ 5 with larger 8-bit samples.
