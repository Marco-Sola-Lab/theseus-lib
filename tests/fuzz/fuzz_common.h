/*
 * Code shared by the fuzzers (fuzz_s2g.cpp and fuzz_msa.cpp).
 *
 * A fuzzer generates random inputs and compares Theseus with a simple, baseline implementation.
 * This header contains everything that does not depend on which aligner is being tested:
 *
 *   1. Environment variables    env_set / env_int
 *   2. Penalties                Pen, make_penalties, gap_cost, describe
 *   3. Random inputs            random_dna, mutate
 *   4. Plain graph              SimpleGraph (a copy of the graph that the reference can read)
 *   5. Checks on an alignment   reference (optimal score), score_ops, replay
 *   6. Reporting                report_failure, print_summary
 *
 * All the aligners are expected to work "end-to-end on the query": the whole query must be
 * aligned, starting at a given (node, offset) of the graph, and it can end anywhere in the graph.
 */
#pragma once

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "theseus/theseus_aligner.h"
#include "theseus/theseus_msa_aligner.h"

using namespace theseus;

// ---------------------------------------------------------------------------------------------
// 1. Environment variables
// ---------------------------------------------------------------------------------------------

// True if the environment variable is defined (e.g. VERB=1, BIG=1...).
inline bool env_set(const char *name) { return std::getenv(name) != nullptr; }

// Integer value of an environment variable, or `fallback` if it is not defined.
inline int env_int(const char *name, int fallback) {
  const char *v = std::getenv(name);
  return v ? std::atoi(v) : fallback;
}

// ---------------------------------------------------------------------------------------------
// 2. Penalties
// ---------------------------------------------------------------------------------------------

// Alignment penalties, as plain numbers (Theseus' own Penalties class is built from them).
// Gap of length k costs  o1 + k*e1,  or  min(o1 + k*e1, o2 + k*e2)  with dual affine gaps.
struct Pen {
  int m, x;     // match, mismatch
  int o1, e1;   // gap open / extend (first piece)
  int o2, e2;   // gap open / extend (second piece, only used if `dual`)
  bool dual;    // dual affine gaps
};

// Cost of a whole gap of length k.
inline long gap_cost(const Pen &p, int k) {
  long c = p.o1 + (long)k * p.e1;
  if (p.dual) c = std::min(c, (long)p.o2 + (long)k * p.e2);
  return c;
}

// Penalties object that Theseus understands.
inline Penalties to_theseus(const Pen &p) {
  return p.dual ? Penalties(p.m, p.x, p.o1, p.e1, p.o2, p.e2) : Penalties(p.m, p.x, p.o1, p.e1);
}

// Chooses the penalties of a trial.
//   Single affine: one of the preset configs, picked randomly (or fixed with CFG=<0..5>).
//   Dual affine:   one of three configs, picked randomly.
// Overrides: PEN="m,x,o,e" replaces the single affine values, DM=<n> sets the match score in dual mode.
// NOTE: the order of the random draws is part of the seed -> do not reorder them.
inline Pen make_penalties(std::mt19937 &rng, bool dual) {
  static const Pen presets[] = {
      {0, 2, 3, 1, 0, 0, false}, {0, 4, 6, 2, 0, 0, false}, {0, 3, 5, 2, 0, 0, false},
      {0, 1, 2, 1, 0, 0, false}, {1, 4, 4, 2, 0, 0, false}, {-1, 2, 3, 1, 0, 0, false}};
  int cfg = env_set("CFG") ? env_int("CFG", 0) : (int)(rng() % 4);
  Pen p = presets[cfg];
  if (env_set("PEN")) std::sscanf(std::getenv("PEN"), "%d,%d,%d,%d", &p.m, &p.x, &p.o1, &p.e1);

  if (dual) {
    p.dual = true;
    int c2 = rng() % 3;
    if (env_set("DM")) p.m = env_int("DM", 0);
    if (c2 == 0)      { p.o1 = 3; p.e1 = 2; p.o2 = 8;  p.e2 = 1; p.x = 3; }
    else if (c2 == 1) { p.o1 = 4; p.e1 = 3; p.o2 = 10; p.e2 = 1; p.x = 4; }
    else              { p.o1 = 2; p.e1 = 2; p.o2 = 6;  p.e2 = 1; p.x = 2; }
  }
  return p;
}

// "(m,x,o1,e1,o2,e2)", used in failure messages.
inline std::string describe(const Pen &p) {
  std::ostringstream os;
  os << "(" << p.m << "," << p.x << "," << p.o1 << "," << p.e1 << "," << p.o2 << "," << p.e2 << ")";
  return os.str();
}

// ---------------------------------------------------------------------------------------------
// 3. Random inputs
// ---------------------------------------------------------------------------------------------

inline char random_base(std::mt19937 &rng) { return "ACGT"[rng() % 4]; }

inline std::string random_dna(std::mt19937 &rng, int len) {
  std::string s;
  for (int i = 0; i < len; ++i) s += random_base(rng);
  return s;
}

// Copy of `s` with random errors. Each character, with probability `rate`, gets one of:
//   - a mismatch (replaced by a random base)            probability rate/3
//   - a deletion of 1..maxgap characters (skipped)      probability rate/3
//   - an insertion of 1..maxgap random bases after it   probability rate/3
// The result is never empty.
inline std::string mutate(std::mt19937 &rng, const std::string &s, double rate, int maxgap) {
  std::string r;
  for (size_t i = 0; i < s.size(); ++i) {
    double u = std::uniform_real_distribution<>(0, 1)(rng);
    if (u < rate / 3) {
      r += random_base(rng);
    } else if (u < 2 * rate / 3) {
      i += rng() % maxgap;
    } else if (u < rate) {
      r += s[i];
      for (int k = 1 + rng() % maxgap; k > 0; --k) r += random_base(rng);
    } else {
      r += s[i];
    }
  }
  if (r.empty()) r = "A";
  return r;
}

// ---------------------------------------------------------------------------------------------
// 4. Plain graph
// ---------------------------------------------------------------------------------------------

// Sequence graph stored with plain vectors: vertex v has the sequence seq[v] and the
// successors out[v]. Fuzzers build it first, and then copy it to a Theseus graph (or read it
// from a GFA), so the reference never depends on Theseus' own graph class.
struct SimpleGraph {
  std::vector<std::string> seq;
  std::vector<std::vector<int>> out;

  void resize(size_t n) { seq.resize(n); out.resize(n); }
  size_t size() const { return seq.size(); }

  // Adds the edge a -> b unless it already exists.
  void add_edge(int a, int b) {
    if (std::find(out[a].begin(), out[a].end(), b) == out[a].end()) out[a].push_back(b);
  }
  bool has_edge(int a, int b) const { return std::find(out[a].begin(), out[a].end(), b) != out[a].end(); }
};

// ---------------------------------------------------------------------------------------------
// 5. Checks on an alignment
// ---------------------------------------------------------------------------------------------

constexpr long INF = 1000000000L;

// Optimal alignment score computed with Dijkstra (all costs are >= 0 except possibly the match
// score, which is why the fuzzers keep it at 0 unless asked otherwise).
//
// Search space: (position in the graph, query offset, state), where
//   - position: every vertex v owns the positions base[v] .. base[v]+len(v); position base[v]+len(v)
//     means "all the vertex is consumed", and it is where the edges to the next vertices start.
//   - state: the last operation (needed for affine gaps: continuing a gap is cheaper than opening one).
// The query is aligned from the start (vertex, offset) and may end at any position of the graph.
// Returns INF if it cannot be aligned.
inline long reference(const SimpleGraph &g, int start, int soff, const std::string &q, const Pen &p) {
  // States: last operation was ...
  enum State { M = 0, I1, D1, I2, D2 };   // match/mismatch, insertion/deletion of gap piece 1 and 2
  constexpr int NSTATES = 5;
  // Insertion = graph character without query character. Deletion = query character without graph character.

  const int n = q.size(), N = g.size();

  // Flat position of each vertex.
  std::vector<int> base(N + 1, 0);
  for (int v = 0; v < N; ++v) base[v + 1] = base[v] + g.seq[v].size() + 1;
  const int P = base[N];

  auto node_id = [&](int pos, int i, int s) { return ((long)pos * (n + 1) + i) * NSTATES + s; };
  std::vector<long> dist((size_t)P * (n + 1) * NSTATES, INF);
  typedef std::pair<long, long> Entry;   // (distance, node id)
  std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

  // Relax: reach (pos, i, s) with cost d.
  auto push = [&](int pos, int i, int s, long d) {
    long k = node_id(pos, i, s);
    if (d < dist[k]) { dist[k] = d; pq.push({d, k}); }
  };
  // Cost of a gap step in state `piece_state` coming from state `from`: extend if it continues the gap, else open it.
  auto step = [&](int from, int piece_state, int o, int e) { return from == piece_state ? e : o + e; };

  push(base[start] + soff, 0, M, 0);
  long best = INF;

  while (!pq.empty()) {
    auto [d, k] = pq.top();
    pq.pop();
    if (d > dist[k]) continue;   // stale queue entry

    // Decode the node.
    int s = k % NSTATES;
    long r = k / NSTATES;
    int i = r % (n + 1);
    int pos = r / (n + 1);
    int v = std::upper_bound(base.begin(), base.end(), pos) - base.begin() - 1;
    int j = pos - base[v];
    int len = g.seq[v].size();

    if (i == n) best = std::min(best, d);   // whole query consumed: valid end (anywhere in the graph)

    // End of the vertex: free move to the start of each successor.
    if (j == len)
      for (int sv : g.out[v]) push(base[sv], i, s, d);

    if (j < len) {
      // Match / mismatch: consumes a query char and a graph char.
      if (i < n) push(pos + 1, i + 1, M, d + (q[i] == g.seq[v][j] ? p.m : p.x));
      // Insertion: consumes only a graph char.
      push(pos + 1, i, I1, d + step(s, I1, p.o1, p.e1));
      if (p.dual) push(pos + 1, i, I2, d + step(s, I2, p.o2, p.e2));
    }
    // Deletion: consumes only a query char.
    if (i < n) {
      push(pos, i + 1, D1, d + step(s, D1, p.o1, p.e1));
      if (p.dual) push(pos, i + 1, D2, d + step(s, D2, p.o2, p.e2));
    }
  }
  return best;
}

// Score of an edit script (string of M/X/I/D) under the penalties `p`.
// Each run of equal gap operations is one gap, so it is charged with gap_cost once.
inline long score_ops(const std::vector<char> &ops, const Pen &p) {
  long s = 0;
  size_t i = 0;
  while (i < ops.size()) {
    char c = ops[i];
    if (c == 'M') { s += p.m; ++i; }
    else if (c == 'X') { s += p.x; ++i; }
    else {
      size_t j = i;
      while (j < ops.size() && ops[j] == c) ++j;
      s += gap_cost(p, j - i);
      i = j;
    }
  }
  return s;
}

// True if consecutive vertices of the path are connected by edges, and it starts at `start`.
inline bool path_is_valid(const SimpleGraph &g, int start, const std::vector<Graph::NodeId> &path, std::string &why) {
  if (path.empty() || (int)path[0] != start) { why = "path does not start at start"; return false; }
  for (size_t k = 0; k + 1 < path.size(); ++k)
    if (!g.has_edge(path[k], path[k + 1])) { why = "invalid path edge"; return false; }
  return true;
}

// Replays the alignment on the graph and the query to check it is real:
//   - the path follows existing edges,
//   - M/X agree with the actual characters (M only if equal, X only if different),
//   - the operations do not run past the end of a vertex or of the query,
//   - the whole query is consumed.
// On failure returns false and explains it in `why`.
inline bool replay(const SimpleGraph &g, int start, int soff, const std::string &q, const Alignment &a, std::string &why) {
  if (!path_is_valid(g, start, a.path, why)) return false;

  size_t pi = 0;       // index in the path of the current vertex
  size_t j = soff;     // offset inside the current vertex
  size_t qi = 0;       // offset inside the query
  auto cur = [&]() -> const std::string & { return g.seq[a.path[pi]]; };
  // When the current vertex is fully consumed, move to the next one in the path.
  auto advance = [&]() {
    while (pi + 1 < a.path.size() && j == cur().size()) { ++pi; j = 0; }
  };

  for (char c : a.edit_op) {
    advance();
    if (c == 'M' || c == 'X') {
      if (qi >= q.size() || j >= cur().size()) { why = "ops overflow"; return false; }
      if ((q[qi] == cur()[j]) != (c == 'M')) { why = "M/X mismatch with chars"; return false; }
      ++qi; ++j;
    } else if (c == 'I') {   // graph char only
      if (j >= cur().size()) { why = "I overflow"; return false; }
      ++j;
    } else {                 // 'D': query char only
      ++qi;
    }
  }
  if (qi != q.size()) { why = "query not consumed"; return false; }
  return true;
}

// ---------------------------------------------------------------------------------------------
// 6. Reporting
// ---------------------------------------------------------------------------------------------

// Prints a failing trial (only the first `max_shown` failures, to keep the output readable).
inline void report_failure(int fails, int max_shown, const std::string &message) {
  if (fails <= max_shown) std::cout << "FAIL " << message << "\n";
}

inline void print_summary(int fails, int trials) { std::cout << "failures: " << fails << " / " << trials << "\n"; }
