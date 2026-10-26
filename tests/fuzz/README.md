# Fuzzers (random test generation)

Unit tests only cover the cases we thought about. The fuzzers generate thousands of random inputs and check Theseus against a simple, obviously-correct implementation (a *differential* test). They are **not** built by CMake, and are meant to be run by hand after changing the aligners.

| File | What it tests |
|------|---------------|
| `fuzz_s2g.cpp` | Sequence-to-graph aligner (`TheseusAligner`). |
| `fuzz_msa.cpp` | MSA aligner (`TheseusMSA`). |
| `fuzz_common.h` | Code shared by both: environment options, penalties, random inputs (`mutate`), the reference aligner, and the score and replay checks. Each fuzzer is organised as input generation → run Theseus → checks (`run_trial`), so a new fuzzer only needs to provide those. |

## Building and running

```
g++ -std=c++23 -O1 -Iinclude tests/fuzz/fuzz_s2g.cpp build/libtheseus.a -o fuzz_s2g
g++ -std=c++23 -O1 -Iinclude tests/fuzz/fuzz_msa.cpp build/libtheseus.a -o fuzz_msa

./fuzz_s2g <trials> <dual 0|1> <cycles 0|1>
./fuzz_msa <trials> <dual 0|1>
```

* `trials`: number of random test cases.
* `dual`: use dual affine gap penalties (two gap pieces) instead of single affine.
* `cycles`: allow cycles in the random graphs (`fuzz_s2g` only).

The program prints the first failing trials and a final `failures: X / N` line. Zero failures means everything agreed.

Optional environment variables: `CFG` (fix the penalty configuration), `PEN="m,x,o,e"` (custom penalties), `BIG` (larger graphs), `LONG` (long vertices and queries), `EMPTY` (empty vertices), `DM` (match score in dual mode), `VERB` (print each trial), `T0` (first trial to run), and for `fuzz_msa`: `TMPF` (temporary GFA file) and `NOPOA` (skip the MSA integrity check). `EMPTY=1` with cycles is expected to hang (cycles of empty vertices).

## What is done

Each trial is generated from a seed derived from the trial number, so **every failure is reproducible**: re-run it with `T0=<trial> VERB=1`.

**`fuzz_s2g`**
1. Generate a random graph: a backbone chain plus random extra edges (and cycles if enabled).
2. Pick a start position, walk the graph to get a path, and mutate it (mismatches, insertions, deletions) to create the query.
3. Align the query with Theseus, and compute the optimal score with the reference.
4. Check that:
   * the score equals the reference score,
   * the edit script (CIGAR) is valid: it follows real graph edges, `M`/`X` match the actual characters, and the whole query is consumed (*replay*),
   * the reported score equals the score recomputed from the CIGAR,
   * the score-only mode returns the same score.

**`fuzz_msa`**
1. Generate a random base sequence and build a POA graph by adding several mutated versions of it.
2. Before adding each sequence, dump the current graph as GFA and align the new sequence against it with the reference.
3. Check the same score and replay conditions as above, plus that the MSA printed by Theseus reproduces every added sequence once the gaps are removed.

## Rationale

* **Reference aligner:** a plain Dijkstra over (graph position, query offset, gap state).
* **Score and replay:** the same score could be reached by a wrong alignment, so the CIGAR is also checked against the graph and the query.
* **Small inputs and several penalty sets:** small graphs and short sequences make failures easy to debug, while many random trials and different penalties (including dual affine) cover corner cases such as self-loops, cycles, long gaps, and vertices of different sizes.
