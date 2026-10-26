/*
 * Differential fuzzer for the MSA aligner (TheseusMSA).
 *
 * Randomized test, NOT built by CMake. Build (from the repository root):
 *   g++ -std=c++23 -O1 -Iinclude tests/fuzz/fuzz_msa.cpp build/libtheseus.a -o fuzz_msa
 *
 * Usage: ./fuzz_msa <trials> <dual 0|1>
 *   env: CFG, PEN="m,x,o,e", DM  (penalties, see make_penalties in fuzz_common.h)
 *        TMPF   temporary GFA file (default /tmp/fuzz_msa.gfa)
 *        NOPOA  skip the check of the printed MSA
 *        VERB   print each trial
 *        T0     first trial to run (to reproduce a failure)
 *
 * Every trial builds a POA graph with TheseusMSA: a random base sequence, then several mutated
 * copies of it. Before adding each copy, the current graph is dumped as GFA and the copy is
 * aligned against it with the reference, starting at the source (vertex 0) and ending anywhere.
 * The checks are those of fuzz_s2g plus the integrity of the printed MSA.
 */
#include "fuzz_common.h"

// Reads the graph (S and L lines) of a GFA file. Vertex names must be integers.
SimpleGraph read_gfa(const std::string &file) {
  std::ifstream in(file);
  std::vector<std::string> seqs;
  std::vector<std::pair<int, int>> edges;

  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::string type;
    ls >> type;
    if (type == "S") {   // S <id> <sequence>
      int id;
      std::string s;
      ls >> id >> s;
      if ((int)seqs.size() <= id) seqs.resize(id + 1);
      seqs[id] = s;
    } else if (type == "L") {   // L <from> <orient> <to> ...
      int a, b;
      std::string orient;
      ls >> a >> orient >> b;
      edges.push_back({a, b});
    }
  }

  SimpleGraph g;
  g.seq = seqs;
  g.out.resize(seqs.size());
  for (auto &e : edges) g.out[e.first].push_back(e.second);
  return g;
}

// The printed MSA (FASTA-like: '>' header lines, then rows with '-' for gaps) must contain
// every sequence that was added, once the gaps are removed. Returns false and sets `why` if not.
bool msa_reproduces_sequences(TheseusMSA &msa, const std::vector<std::string> &added, std::string &why) {
  std::ostringstream os;
  msa.print_as_msa(os);
  std::istringstream is(os.str());

  // Rows without gaps, in order.
  std::vector<std::string> rows;
  std::string line;
  while (std::getline(is, line)) {
    if (!line.empty() && line[0] == '>') rows.push_back("");
    else if (!rows.empty())
      for (char c : line)
        if (c != '-') rows.back() += c;
  }

  // One row per added sequence, plus one extra row (printed after them).
  if (rows.size() != added.size() + 1) {
    why = "MSA rows=" + std::to_string(rows.size());
    return false;
  }
  for (size_t r = 0; r < added.size(); ++r)
    if (rows[r] != added[r]) {
      why = "MSA row " + std::to_string(r) + " differs from its sequence";
      return false;
    }
  return true;
}

// Runs a trial. Returns true if everything agrees; otherwise `message` describes the failure.
bool run_trial(int t, bool dual, const std::string &tmp_gfa, std::string &message) {
  std::mt19937 rng(t * 7919 + 13);
  Pen p = make_penalties(rng, dual);

  // The base sequence is the first one of the POA graph.
  std::string base = random_dna(rng, 20 + rng() % 100);
  Heuristics h;
  TheseusMSA msa(to_theseus(p), h, base, 1);
  std::vector<std::string> added = {base};

  int nseq = 2 + rng() % 5;
  for (int s = 0; s < nseq; ++s) {
    std::string q = mutate(rng, base, 0.05 + (rng() % 20) / 100.0, 1 + rng() % (dual ? 20 : 5));

    // Reference: optimal score against the graph as it is now (before adding q).
    { std::ofstream o(tmp_gfa); msa.print_as_gfa(o); }
    SimpleGraph g = read_gfa(tmp_gfa);
    long ref = reference(g, 0, 0, q, p);

    // Theseus: align q and add it to the graph.
    Alignment a = msa.align(q);
    std::string why;
    long got = score_ops(a.edit_op, p);
    bool ok = replay(g, 0, 0, q, a, why);
    added.push_back(q);

    if (ok && !env_set("NOPOA")) ok = msa_reproduces_sequences(msa, added, why);

    if (!ok || got != ref) {
      message = "trial " + std::to_string(t) + " seq " + std::to_string(s) + " cfg" + describe(p) +
                " ref=" + std::to_string(ref) + " got=" + std::to_string(got) + " " + why;
      return false;   // the graph is already wrong: following sequences would not tell us more
    }
  }
  return true;
}

int main(int argc, char **argv) {
  int trials = argc > 1 ? atoi(argv[1]) : 200;
  bool dual = argc > 2 ? atoi(argv[2]) : false;
  std::string tmp_gfa = env_set("TMPF") ? std::getenv("TMPF") : "/tmp/fuzz_msa.gfa";

  int fails = 0;
  for (int t = env_int("T0", 0); t < trials; ++t) {
    if (env_set("VERB")) std::cerr << "t=" << t << std::endl;
    std::string message;
    if (!run_trial(t, dual, tmp_gfa, message)) report_failure(++fails, 6, message);
  }
  print_summary(fails, trials);
}
