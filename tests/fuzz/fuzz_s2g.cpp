/*
 * Differential fuzzer for the sequence-to-graph aligner (TheseusAligner).
 *
 * Randomized test, NOT built by CMake. Build (from the repository root):
 *   g++ -std=c++23 -O1 -Iinclude tests/fuzz/fuzz_s2g.cpp build/libtheseus.a -o fuzz_s2g
 *
 * Usage: ./fuzz_s2g <trials> <dual 0|1> <cycles 0|1>
 *   env: CFG, PEN="m,x,o,e", DM  (penalties, see make_penalties in fuzz_common.h)
 *        BIG    larger graphs
 *        LONG   some vertices up to 300 bases, and longer queries
 *        EMPTY  vertices can be empty (hangs with cycles=1: cycles of empty vertices)
 *        VERB   print each trial (graph, query...)
 *        T0     first trial to run (to reproduce a failure)
 *
 * Every trial:
 *   1. builds a random graph (backbone chain + random extra edges, cycles if enabled),
 *   2. takes a random walk on it and mutates it to get the query,
 *   3. aligns the query end-to-end from (start vertex, offset), free end in the graph,
 *   4. compares with the reference (Dijkstra) and replays the CIGAR, see fuzz_common.h.
 */
#include "fuzz_common.h"

// Options read from the environment (constant for the whole run).
struct Options {
  bool big = env_set("BIG");
  bool long_vertices = env_set("LONG");
  bool empty_vertices = env_set("EMPTY");
  bool verbose = env_set("VERB");
};

// ---------------------------------------------------------------------------------------------
// Input generation
// NOTE: the order of the random draws defines what a trial contains (seeds are reproducible),
// so keep it when editing.
// ---------------------------------------------------------------------------------------------

// Random graph: vertices with random sequences, a backbone chain v -> v+1 so everything is
// reachable, and some random extra edges. Edges going backwards (cycles, self-loops) are only
// added if `cycles` is set.
SimpleGraph random_graph(std::mt19937 &rng, bool cycles, const Options &opt) {
  int N = 2 + rng() % (opt.big ? 30 : 9);
  SimpleGraph g;
  g.resize(N);

  for (int v = 0; v < N; ++v) {
    // Length: at least 1 (unless EMPTY), up to 7. With LONG, one vertex out of three is up to 300 bases.
    int len = (opt.empty_vertices ? 0 : 1) + (opt.long_vertices && rng() % 3 == 0 ? rng() % 300 : rng() % 8);
    g.seq[v] = random_dna(rng, len);
  }

  for (int v = 0; v + 1 < N; ++v) g.add_edge(v, v + 1);
  int extra = rng() % (N + 2);
  for (int e = 0; e < extra; ++e) {
    int a = rng() % N, b = rng() % N;
    if (cycles || b > a) g.add_edge(a, b);
  }
  return g;
}

// Copy of the graph in Theseus' own graph class.
Graph to_theseus_graph(const SimpleGraph &g) {
  Graph gr;
  for (const std::string &s : g.seq) (void)gr.add_node(s);
  for (size_t v = 0; v < g.size(); ++v)
    for (int w : g.out[v]) gr.add_edge(v, w);
  return gr;
}

// Spells the sequence of a random path that starts at (start, soff), until it has at least
// `target` characters or the path ends (a vertex without successors, or 1000 steps).
std::string random_walk(std::mt19937 &rng, const SimpleGraph &g, int start, int soff, int target) {
  std::string walk;
  int v = start, j = soff;
  for (int steps = 0; (int)walk.size() < target && steps < 1000; ++steps) {
    if (j == (int)g.seq[v].size()) {   // end of vertex: go to a random successor
      if (g.out[v].empty()) break;
      v = g.out[v][rng() % g.out[v].size()];
      j = 0;
    } else {
      walk += g.seq[v][j++];
    }
  }
  return walk;
}

// ---------------------------------------------------------------------------------------------
// Checks
// ---------------------------------------------------------------------------------------------

void print_graph(const SimpleGraph &g) {
  for (size_t v = 0; v < g.size(); ++v) {
    std::cerr << v << ":\"" << g.seq[v] << "\" ->";
    for (int w : g.out[v]) std::cerr << " " << w;
    std::cerr << "\n";
  }
}

// Runs a trial. Returns true if everything agrees; otherwise `message` describes the failure.
bool run_trial(int t, bool dual, bool cycles, const Options &opt, std::string &message) {
  std::mt19937 rng(t * 104729 + 7);
  Pen p = make_penalties(rng, dual);

  // Random graph, start position and query.
  SimpleGraph g = random_graph(rng, cycles, opt);
  int N = g.size();
  int start = rng() % N;
  int soff = rng() % std::max<size_t>(1, g.seq[start].size());
  int target = 5 + rng() % (opt.long_vertices ? 400 : 40);
  std::string walk = random_walk(rng, g, start, soff, target);
  double rate = 0.05 + (rng() % 25) / 100.0;
  int maxgap = dual ? 14 : 4;
  std::string q = mutate(rng, walk, rate, maxgap);

  if (opt.verbose) print_graph(g);
  long ref = reference(g, start, soff, q, p);
  if (opt.verbose)
    std::cerr << "ref done " << ref << " q=" << q << " start=" << start << "," << soff << " N=" << N << std::endl;

  // Align with Theseus (full mode: score + CIGAR + path).
  Heuristics h;
  TheseusAligner aligner(to_theseus(p), h, to_theseus_graph(g));
  Graph::NodeId start_node = start;
  Alignment a = aligner.align(q, start_node, soff);

  std::string why;
  bool ok = replay(g, start, soff, q, a, why);
  long got = score_ops(a.edit_op, p);
  if (ok && a.score != got) {   // the reported score must match the CIGAR
    ok = false;
    why = "reported score " + std::to_string(a.score) + " != CIGAR score";
  }

  // Score-only mode must give the same score as the reference.
  aligner.set_alignment_mode(AlignmentMode::ScoreOnly);
  Alignment score_only = aligner.align(q, start_node, soff);
  if (ok && score_only.score != ref) {
    ok = false;
    why = "score-only " + std::to_string(score_only.score);
  }

  if (ok && got == ref) return true;
  message = "trial " + std::to_string(t) + " cfg" + describe(p) + " ref=" + std::to_string(ref) +
            " got=" + std::to_string(got) + " " + why + " qlen=" + std::to_string(q.size()) +
            " N=" + std::to_string(N);
  return false;
}

int main(int argc, char **argv) {
  int trials = argc > 1 ? atoi(argv[1]) : 200;
  bool dual = argc > 2 ? atoi(argv[2]) : false;
  bool cycles = argc > 3 ? atoi(argv[3]) : true;
  Options opt;

  int fails = 0;
  for (int t = env_int("T0", 0); t < trials; ++t) {
    if (opt.verbose) std::cerr << "t=" << t << std::endl;
    std::string message;
    if (!run_trial(t, dual, cycles, opt, message)) report_failure(++fails, 8, message);
  }
  print_summary(fails, trials);
}
