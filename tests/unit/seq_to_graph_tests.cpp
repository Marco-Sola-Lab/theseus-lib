/*
 *                             The MIT License
 *
 * Copyright (c) 2024 by Albert Jimenez-Blanco
 *
 * This file is part of #################### Theseus Library ####################.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */


#include "../doctest.h"

#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "../../include/theseus/graph.h"
#include "../../include/theseus/alignment.h"
#include "../../include/theseus/penalties.h"
#include "../../include/theseus/theseus_aligner.h"


using NodeId = theseus::Graph::NodeId;

// Score of a CIGAR under the user penalties (a gap is a maximal run of I or D,
// and with dual-affine gaps it costs the cheapest of the two pieces)
static int cigar_score(const std::vector<char> &ops, const theseus::Penalties &p) {
    int score = 0;
    size_t i = 0;
    while (i < ops.size()) {
        if (ops[i] == 'M') { score += p.match(); ++i; }
        else if (ops[i] == 'X') { score += p.mism(); ++i; }
        else {
            size_t j = i;
            while (j < ops.size() && ops[j] == ops[i]) ++j;
            int len = j - i;
            int gap = p.gapo() + len * p.gape();
            if (p.type() == theseus::Penalties::Type::DualAffine) {
                gap = std::min(gap, p.gapo2() + len * p.gape2());
            }
            score += gap;
            i = j;
        }
    }
    return score;
}

TEST_CASE("Check sequence-to-graph aligner") {
    SUBCASE("Correct alignment of sequences against a graph with a cycle") {
        // Reference graph
        theseus::Graph G;
        // Add nodes
        NodeId n1 = G.add_node("ACTTAG");
        NodeId n2 = G.add_node("ACA");
        NodeId n3 = G.add_node("T");
        NodeId n4 = G.add_node("GTACTT");
        // Add edges
        G.add_edge(n1, n2);
        G.add_edge(n1, n3);
        G.add_edge(n2, n4);
        G.add_edge(n3, n4);
        G.add_edge(n4, n1);

        // Sequences for alignment
        std::vector<std::string> sequences = {
            "TAGACAGTACT",   // Perfect match
            "TAGACAGGACT",   // One mismatch
            "ACAGTACTTACT",  // Perfect match with a cycle
            "AACAGTACTTACT", // A deletion with a cycle
            "ACAGTATTACT"    // An insertion with a cycle
        };

        // Starting offsets and vertices
        std::vector<int> start_offsets = {3, 3, 0, 0, 0};
        std::vector<NodeId> start_vertices = {n1, n1, n2, n2, n2};

        // Expected CIGARs
        std::vector<std::vector<char>> expected_cigars = {
            {'M','M','M','M','M','M','M','M','M','M','M'},
            {'M','M','M','M','M','M','M','X','M','M','M'},
            {'M','M','M','M','M','M','M','M','M','M','M','M'},
            {'M','D','M','M','M','M','M','M','M','M','M','M','M'},
            {'M','M','M','M','M','M','I','M','M','M','M','M'}
        };

        // Expected paths
        std::vector<std::vector<NodeId>> expected_paths = {
            {n1, n2, n4},
            {n1, n2, n4},
            {n2, n4, n1},
            {n2, n4, n1},
            {n2, n4, n1}
        };

        // Expected scores
        std::vector<int> expected_scores = {0, 2, 0, 4, 4};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                            // Create penalties object
        theseus::Heuristics heuristics;                        // Create heuristics object
        theseus::TheseusAligner aligner(penalties, heuristics, std::move(G));       // Create aligner

        // Align sequences and check results
        for (int i = 0; i < start_offsets.size(); ++i) {
            theseus::Alignment alignment = aligner.align(
                sequences[i], start_vertices[i], start_offsets[i], false, false
            );

            CHECK(alignment.score == expected_scores[i]); // Check score
            CHECK(alignment.edit_op == expected_cigars[i]);        // Check CIGAR
            CHECK(alignment.path == expected_paths[i]); // Check path
        }
    }
}

// Regression tests. Expected scores computed with a reference Dijkstra/Gotoh over
// the graph (query aligned end to end, free end in the graph).
TEST_CASE("Check sequence-to-graph aligner regressions") {
    SUBCASE("Insertion crossing a vertex that jumps to itself (self-loop)") {
        theseus::Graph G;
        NodeId n1 = G.add_node("T");
        NodeId n2 = G.add_node("TAGATGAG");
        G.add_edge(n1, n1);
        G.add_edge(n1, n2);

        theseus::Penalties penalties(0, 1, 2, 1);
        theseus::Heuristics heuristics;
        theseus::TheseusAligner aligner(penalties, heuristics, std::move(G));

        NodeId start = n1;
        theseus::Alignment alignment = aligner.align("TTTAATAGATGAGCGAG", start, 0);
        CHECK(alignment.score == 8);
    }

    SUBCASE("Match score different from 0") {
        theseus::Graph G;
        std::vector<NodeId> n;
        for (std::string s : {"TTA", "TCCG", "AAGAACG", "GCCACCG", "T", "TCGA", "AT"}) n.push_back(G.add_node(s));
        G.add_edge(n[0], n[1]); G.add_edge(n[1], n[2]); G.add_edge(n[1], n[6]);
        G.add_edge(n[2], n[3]); G.add_edge(n[3], n[4]); G.add_edge(n[4], n[5]); G.add_edge(n[5], n[6]);

        theseus::Penalties penalties(1, 4, 4, 2);
        theseus::Heuristics heuristics;
        theseus::TheseusAligner aligner(penalties, heuristics, std::move(G));

        NodeId start = n[0];
        theseus::Alignment alignment = aligner.align("AAACCGAAAAACCGGCCTGTTCG", start, 1);
        CHECK(alignment.score == 46);
    }

    SUBCASE("Heuristic flags do not change the direction of the alignment") {
        theseus::Graph G;
        NodeId n1 = G.add_node("ACGTTGCA");
        NodeId n2 = G.add_node("AGCTTAGG");
        G.add_edge(n1, n2);
        theseus::Penalties penalties(0, 2, 3, 1);
        theseus::Heuristics heuristics;
        theseus::TheseusAligner aligner(penalties, heuristics, std::move(G));

        NodeId start = n1;
        // The flags are forwarded as density drop / lag pruning (they used to be
        // received as "ends free" and "reverse alignment", which walks the graph
        // backwards and cannot reach the second vertex)
        theseus::Alignment alignment = aligner.align("ACGTTGCAAGCTTAGG", start, 0, false, true);
        CHECK(alignment.score == 0);
        CHECK(alignment.path == std::vector<NodeId>({n1, n2}));
        alignment = aligner.align("ACGTTGCAAGCTTAGG", start, 0, true, false);
        CHECK(alignment.score == 0);
        CHECK(alignment.path == std::vector<NodeId>({n1, n2}));
    }
}


TEST_CASE("Check sequence-to-graph aligner with dual affine gaps") {
    // Gap cost = min(6 + 2k, 12 + k)
    std::string seq = "ACGTTGCAAGCTTAGGCATCGATCCGATTACGGATCAGTC";
    theseus::Penalties penalties(0, 4, 6, 2, 12, 1);
    theseus::Heuristics heuristics;

    SUBCASE("Penalties are validated") {
        theseus::Graph G;
        (void)G.add_node(seq);
        // The second piece must have a higher open and a lower extension penalty
        CHECK_THROWS(theseus::TheseusAligner(theseus::Penalties(0, 4, 6, 2, 5, 1), heuristics, G));
        CHECK_THROWS(theseus::TheseusAligner(theseus::Penalties(0, 4, 6, 2, 12, 3), heuristics, G));
    }

    SUBCASE("Long gaps use the second piece") {
        theseus::Graph G;
        NodeId n1 = G.add_node(seq);
        theseus::TheseusAligner aligner(penalties, heuristics, std::move(G));

        NodeId start = n1;
        // 12 chars missing in the query: 12 + 12 = 24 (instead of 6 + 24 = 30)
        theseus::Alignment alignment = aligner.align(seq.substr(0, 10) + seq.substr(22), start, 0);
        CHECK(alignment.score == 24);
        // 12 extra chars in the query
        alignment = aligner.align(seq.substr(0, 25) + "TTTTTTTTTTTT" + seq.substr(25), start, 0);
        CHECK(alignment.score == 24);
        // Short gaps keep using the first piece
        alignment = aligner.align(seq.substr(0, 10) + seq.substr(13), start, 0);
        CHECK(alignment.score == 12);
    }
}


TEST_CASE("Check the score-only mode") {
    // Graph with a bubble and a cycle
    theseus::Graph G;
    NodeId n1 = G.add_node("ACGTTGCAAGCT");
    NodeId n2 = G.add_node("TAGG");
    NodeId n3 = G.add_node("CATCGA");
    NodeId n4 = G.add_node("TCCGATTACGGATCAGTC");
    G.add_edge(n1, n2);
    G.add_edge(n1, n3);
    G.add_edge(n2, n4);
    G.add_edge(n3, n4);
    G.add_edge(n4, n1);

    std::vector<std::string> sequences = {
        "ACGTTGCAAGCTTAGGTCCGATTACGGATCAGTC",          // Perfect match
        "ACGTTGCTAGCTCATCGATCCGATACGGATCAGTCACGTTGC",  // Errors and a cycle
        "ACGTTGCAAGCTTTTTTTTTTAGGTCCGATTACGG"          // Long insertion
    };

    for (const auto &penalties : {theseus::Penalties(0, 4, 6, 2),
                                  theseus::Penalties(1, 4, 6, 2),
                                  theseus::Penalties(0, 4, 6, 2, 12, 1)}) {
        theseus::Heuristics heuristics;
        theseus::Graph G_copy = G;
        theseus::TheseusAligner aligner(penalties, heuristics, std::move(G_copy));
        for (const auto &seq : sequences) {
            NodeId start = n1;
            theseus::Alignment full = aligner.align(seq, start, 0);
            aligner.set_alignment_mode(theseus::AlignmentMode::ScoreOnly);
            theseus::Alignment score_only = aligner.align(seq, start, 0);
            aligner.set_alignment_mode(theseus::AlignmentMode::Full);

            // Same score (also equal to the score of the CIGAR) and end position
            CHECK(full.score == cigar_score(full.edit_op, penalties));
            CHECK(score_only.score == full.score);
            CHECK(score_only.edit_op.empty());
            CHECK(score_only.path == std::vector<NodeId>({full.path.back()}));
            CHECK(score_only.end_offset == full.end_offset);
        }
    }
}
