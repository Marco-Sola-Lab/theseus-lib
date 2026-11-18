# Tests

This folder contains the tests of Theseus-lib:

* [unit/](unit/): deterministic unit tests, built and run with CMake/CTest.
* [fuzz/](fuzz/): randomized differential tests (see [fuzz/README.md](fuzz/README.md)).

## Unit tests

The unit tests use [doctest](https://github.com/doctest/doctest) (`doctest.h`, with the shared `main` in `doctest_main.cpp`). Each file in `unit/` is compiled into its own executable and registered in CTest, so adding a new test only requires dropping a new `.cpp` file in `unit/`.

| File | What it checks |
|------|----------------|
| `msa_tests.cpp` | MSA aligner: matches, mismatches, insertions and deletions, diverging sequences, several sequences, dual affine gaps, and regressions (weighted majority voting, the MSA reproduces every input sequence). It also checks the per-sequence options: the three end conditions (scores, and whether the sequence is linked to the sink of the POA graph), reverse alignments (linked to the source, mixed with forward ones) and the density drop (partial alignment added to the POA graph, status and score). |
| `seq_to_graph_tests.cpp` | Sequence-to-graph aligner: graphs with cycles, self-loops, non-zero match scores, heuristic flags, penalty validation, dual affine gaps and the score-only mode. |
| `vector.cpp` | Internal vector class: constructors, assignments, resizing, iterators, allocation policies. |

### Running them

They are built by default (`ENABLE_TESTS=ON`). From the build folder:

```
cmake .. -DCMAKE_BUILD_TYPE=Release
make
ctest --output-on-failure
```

The executables are placed in `build/tests/unit/`, and can be run directly (e.g. `./tests/unit/msa_tests -sc="Check MSA aligner"` to run a single test case).

### Rationale

* Every `SUBCASE` is a small, hand-written example where the expected CIGAR and score are known, so a failure points directly to the broken behaviour.
* Each bug that has been fixed gets a *regression* subcase, so it cannot silently come back.
* Unit tests only cover cases we thought about. The fuzzers in [fuzz/](fuzz/) cover the rest.
