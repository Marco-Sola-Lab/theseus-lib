/*
 * Differential fuzzer for the MSA aligner (TheseusMSA).
 *
 * Randomized test, NOT built by CMake. Build (from the repository root):
 *   g++ -std=c++23 -O1 -Iinclude tests/fuzz/fuzz_msa.cpp build/libtheseus.a -o fuzz_msa
 *
 * Usage: ./fuzz_msa <trials> <dual 0|1>
 *   env: CFG, PEN="m,x,o,e", DM  (penalties, see make_penalties in fuzz_common.h)
 *        EC     end condition of the alignments: 0 full query (default), 1 end to end, 2 ends free
 *        REV    align each sequence forward or in reverse (randomly)
 *        DROP   use the density drop heuristic, with sequences that make it drop the alignment
 *        TMPF   temporary GFA file (default /tmp/fuzz_msa.gfa)
 *        NOPOA  skip the check of the printed MSA
 *        VERB   print each trial
 *        T0     first trial to run (to reproduce a failure)
 *
 * Every trial builds a POA graph with TheseusMSA: a random base sequence, then several mutated
 * copies of it. Before adding each copy, the current graph is dumped as GFA and the copy is
 * aligned against it with the reference, starting at the source (vertex 0), or at the sink in a
 * reverse alignment, and ending as the end condition says. With EC=2 or REV some copies are
 * cut, so they only cover the first (forward) or last (reverse) part of the graph.
 * The checks are those of fuzz_s2g, plus the end condition and the integrity of the printed MSA.
 * A dropped alignment (DROP) is partial: only that the partial alignment is valid, its score
 * and that its part of the sequence is in the MSA are checked.
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

// Number of alignments dropped by the density drop heuristic (shown with DROP, to know that
// the partial alignments are really tested).
int dropped_alignments = 0;

// Part of the sequence q that is in the MSA when the alignment consumes `consumed` characters:
// the first ones in a forward alignment, the last ones in a reverse alignment.
std::string aligned_part(const std::string &q, size_t consumed, bool reversed) {
  return reversed ? q.substr(q.size() - consumed) : q.substr(0, consumed);
}

// Sequence to add to the POA graph: a mutated copy of the base sequence. With the options of
// the fuzzer it can be cut (it only covers the beginning of the graph, or the end if it is
// aligned in reverse) or have random characters (junk) after its end (before its start if it
// is aligned in reverse), so the density drop heuristic drops its alignment.
std::string make_query(std::mt19937 &rng, const std::string &base, bool dual, bool reversed, int end_cond, bool use_reverse, bool use_drop) {
  std::string q = mutate(rng, base, 0.05 + (rng() % 20) / 100.0, 1 + rng() % (dual ? 20 : 5));
  if ((use_reverse || end_cond != 0) && end_cond != 1 && rng() % 2) {
    size_t len = std::max<size_t>(1, q.size() * (30 + rng() % 61) / 100);
    q = reversed ? q.substr(q.size() - len) : q.substr(0, len);
  }
  if (use_drop && rng() % 2) {
    std::string junk = random_dna(rng, 150 + rng() % 100);
    q = reversed ? junk + q : q + junk;
  }
  return q;
}

// Runs a trial. Returns true if everything agrees; otherwise `message` describes the failure.
bool run_trial(int t, bool dual, const std::string &tmp_gfa, std::string &message) {
  std::mt19937 rng(t * 7919 + 13);
  Pen p = make_penalties(rng, dual);
  const int end_cond = env_int("EC", 0);
  const bool use_reverse = env_set("REV");
  const bool use_drop = env_set("DROP");
  // Vertex of the POA graph where the alignments end if they reach the end of the graph
  // (the sink, or the source in a reverse alignment), and the one where they start
  constexpr int SOURCE = 0, SINK = 2;

  // The base sequence is the first one of the POA graph.
  std::string base = random_dna(rng, (use_drop ? 200 : 20) + rng() % 100);
  Heuristics h;
  TheseusMSA msa(to_theseus(p), h, base, 1);
  std::vector<std::string> added = {base};

  int nseq = 2 + rng() % 5;
  for (int s = 0; s < nseq; ++s) {
    const bool reversed = use_reverse && (rng() % 2);
    std::string q = make_query(rng, base, dual, reversed, end_cond, use_reverse, use_drop);

    // Reference: optimal score against the graph as it is now (before adding q). A reverse
    // alignment is a forward alignment of the reversed query to the reversed graph.
    { std::ofstream o(tmp_gfa); msa.print_as_gfa(o); }
    SimpleGraph g = read_gfa(tmp_gfa);
    SimpleGraph rg = reverse_graph(g);
    std::string rq(q.rbegin(), q.rend());
    const SimpleGraph &walk_graph = reversed ? rg : g;
    const std::string &walk_query = reversed ? rq : q;
    const int walk_start = reversed ? SINK : SOURCE;
    const int walk_end   = reversed ? SOURCE : SINK;
    long ref = reference(walk_graph, walk_start, 0, walk_query, p, end_cond, walk_end);

    // Theseus: align q and add it to the graph.
    MSAOptions options;
    options.reverse       = reversed;
    options.end_condition = (EndCondition)end_cond;
    options.density_drop  = use_drop;
    Alignment a = msa.align(q, 1, options);
    const bool dropped = (a.theseus_status == THESEUS_STATUS_END_UNREACHABLE);
    if (dropped) ++dropped_alignments;

    // The alignment in the direction in which it was computed (a reverse alignment lists the
    // path and the operations from the start of the graph, so it is reversed)
    Alignment walked = a;
    if (reversed) {
      std::reverse(walked.edit_op.begin(), walked.edit_op.end());
      std::reverse(walked.path.begin(), walked.path.end());
    }
    std::string why;
    size_t consumed = 0;
    bool ok = replay(walk_graph, walk_start, 0, walk_query, walked, why, &consumed);
    long got = score_ops(a.edit_op, p);
    if (ok && a.score != got) { ok = false; why = "score of the alignment differs from the one of its CIGAR"; }

    if (ok && !dropped) {
      // Complete alignment: it must end as the end condition says, and be optimal
      const bool query_done = (consumed == q.size());
      const bool end_reached = (walked.path.back() == walk_end);
      if (a.theseus_status != THESEUS_STATUS_ALG_COMPLETED) { ok = false; why = "status is not completed"; }
      else if (end_cond == 0 && !query_done) { ok = false; why = "query not consumed"; }
      else if (end_cond == 1 && !(query_done && end_reached)) { ok = false; why = "end to end not reached"; }
      else if (end_cond == 2 && !(query_done || end_reached)) { ok = false; why = "ends free: no end reached"; }
    }
    else if (ok && !use_drop) { ok = false; why = "alignment dropped without density drop"; }
    added.push_back(aligned_part(q, consumed, reversed));

    if (ok && !env_set("NOPOA")) ok = msa_reproduces_sequences(msa, added, why);

    if (!ok || (!dropped && got != ref)) {
      message = "trial " + std::to_string(t) + " seq " + std::to_string(s) + " cfg" + describe(p) +
                " reverse=" + std::to_string(reversed) + " ref=" + std::to_string(ref) + " got=" + std::to_string(got) + " " + why;
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
  if (env_set("DROP")) std::cout << "dropped alignments: " << dropped_alignments << "\n";
}
