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

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "../../include/theseus/alignment.h"
#include "../../include/theseus/penalties.h"
#include "../../include/theseus/theseus_msa_aligner.h"
#include "../../include/theseus/heuristics.h"
#include "../../include/theseus/graph.h"

using NodeId = theseus::Graph::NodeId;


TEST_CASE("Check MSA aligner") {
    SUBCASE("Correct MSA with a matching sequence") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCGTAAAAGGG";
        std::vector<char> expected_cigar = {'M','M','M','M','M','M','M','M','M','M','M','M','M'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, 1, false);

        // Check if alignment was successful
        CHECK(alignment.score == 0); // Check score
        CHECK(alignment.edit_op == expected_cigar);                  // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));     // Check path
    }

    SUBCASE("Correct MSA with a sequence with a mismatch") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCGTCAAAGGG";
        std::vector<char> expected_cigar = {'M','M','M','M','M','M','X','M','M','M','M','M','M'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, 1, false);

        // Check if alignment was successful
        CHECK(alignment.score == 2); // Check score
        CHECK(alignment.edit_op == expected_cigar);                  // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));        // Check path
    }

    SUBCASE("Correct MSA with an deletion at the end") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCGTAAAAGGGAAA";
        std::vector<char> expected_cigar = {'M','M','M','M','M','M','M','M','M','M','M','M','M', 'D','D','D'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, 1, false);

        // Check if alignment was successful
        CHECK(alignment.score == 6);    // Check score
        CHECK(alignment.edit_op == expected_cigar);                      // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));         // Check path
    }

    SUBCASE("Correct MSA with an deletion at the beginning") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "CATACCCGTAAAAGGG";
        std::vector<char> expected_cigar = {'D','D','D','M','M','M','M','M','M','M','M','M','M','M','M','M'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, 1, false);

        // Check if alignment was successful
        CHECK(alignment.score == 6); // Check score
        CHECK(alignment.edit_op == expected_cigar);                  // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));        // Check path
    }

    SUBCASE("Correct MSA with an insertion in the middle") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCGAAGGG";
        std::vector<char> expected_cigar = {'M','M','M','M','M','I','I','I','M','M','M','M','M'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq);

        // Check if alignment was successful
        CHECK(alignment.score == 6); // Check score
        CHECK(alignment.edit_op == expected_cigar);                  // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));        // Check path
    }

    SUBCASE("Correct MSA with diverging sequence") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCCCATAAGAGGG";

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, 1, false);

        // Check if alignment was successful
        CHECK(alignment.score == 9); // Check score
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));        // Check path
    }

    SUBCASE("Correct MSA with several sequences") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string seq_1       = "ACCCGTCAAAGGG";
        std::string seq_2       = "ACCCGAAGGG";
        std::string seq_3       = "ACCCGTCAAAGGG";
        std::string seq_4       = "ACCCCCATAAGAGGG";

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1); // Create aligner

        // Align and check sequence 1
        theseus::Alignment alignment = aligner.align(seq_1);
        CHECK(alignment.score == 2); // Check score

        // Align and check sequence 2
        alignment = aligner.align(seq_2);
        CHECK(alignment.score == 6); // Check score

        // Align and check sequence 3
        alignment = aligner.align(seq_3);
        CHECK(alignment.score == 0); // Check score

        // Align and check sequence 4
        alignment = aligner.align(seq_4);
        CHECK(alignment.score == 9); // Check score
    }

    SUBCASE("Correct reversed MSA") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string seq_1       = "ACCCGTCAAAGGG";
        std::string seq_2       = "ACCCGAAGGG";
        std::string seq_3       = "GTAAAAGGG";
        // std::string seq_3       = "ACCCGTCAAAGGG";
        // std::string seq_4       = "ACCCCCATAAGAGGG";

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1); // Create aligner

        // Align and check sequence 1
        theseus::Alignment alignment = aligner.align(seq_1, 1, true);
        CHECK(alignment.score == 2); // Check score

        // Align and check sequence 2
        alignment = aligner.align(seq_2, 1, false);
        CHECK(alignment.score == 6); // Check score

        // // Align and check sequence 3
        alignment = aligner.align(seq_3, 1, true);
        CHECK(alignment.score == 7); // Check score

        // // Align and check sequence 4
        // alignment = aligner.align(seq_4);
        // CHECK(alignment.score == 9); // Check score
    }
}

TEST_CASE("Check MSA aligner regressions") {
    SUBCASE("Weighted majority voting counts the last added sequence") {
        // The insertion appears in 1 of 2 sequences: it is not a majority (exactly
        // half), so it must not be part of the consensus
        theseus::Penalties penalties(0, 2, 3, 1);
        theseus::Heuristics heuristics;
        theseus::TheseusMSA aligner(penalties, heuristics, "ACGTACGTACGT", 1);
        aligner.align("ACGTACTTTTGTACGT", 1);

        std::vector<int> weights;
        std::string consensus, gapped_consensus;
        aligner.majority_voting_consensus(weights, consensus, gapped_consensus);
        CHECK(consensus == "ACGTACGTACGT");
    }

    SUBCASE("The MSA reproduces every sequence") {
        theseus::Penalties penalties(0, 2, 3, 1);
        theseus::Heuristics heuristics;
        std::vector<std::string> sequences = {"ACCCGTAAAAGGG", "ACCCGTCAAAGGG", "ACCCGAAGGG", "ACCCCCATAAGAGGG"};
        theseus::TheseusMSA aligner(penalties, heuristics, sequences[0], 1);
        for (size_t i = 1; i < sequences.size(); ++i) aligner.align(sequences[i]);

        std::ostringstream msa;
        aligner.print_as_msa(msa);
        std::istringstream in(msa.str());
        std::string line;
        std::vector<std::string> rows;
        while (std::getline(in, line)) {
            if (line[0] == '>') rows.push_back("");
            else for (char c : line) if (c != '-') rows.back() += c;
        }
        REQUIRE(rows.size() == sequences.size() + 1);   // + consensus
        for (size_t i = 0; i < sequences.size(); ++i) CHECK(rows[i] == sequences[i]);
    }

    SUBCASE("Dual affine gaps") {
        std::string initial_seq = "ACGTTGCAAGCTTAGGCATCGATCCGATTACGGATCAGTC";
        theseus::Penalties penalties(0, 4, 6, 2, 12, 1);
        theseus::Heuristics heuristics;
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);

        theseus::Alignment alignment = aligner.align(initial_seq.substr(0, 10) + initial_seq.substr(22));
        CHECK(alignment.score == 24);
    }
}


// ------------------------- End conditions and options -------------------------

// Rows of the printed MSA without gaps (the last row is the consensus)
static std::vector<std::string> ungapped_rows(theseus::TheseusMSA &aligner) {
    std::ostringstream msa;
    aligner.print_as_msa(msa);
    std::istringstream in(msa.str());
    std::string line;
    std::vector<std::string> rows;
    while (std::getline(in, line)) {
        if (line[0] == '>') rows.push_back("");
        else for (char c : line) if (c != '-') rows.back() += c;
    }
    return rows;
}

// Number of edges of the POA graph (in GFA) that arrive to the sink (vertex 2) or that
// leave the source (vertex 0)
static int count_links(theseus::TheseusMSA &aligner, bool to_sink) {
    std::ostringstream gfa;
    aligner.print_as_gfa(gfa);
    std::istringstream in(gfa.str());
    std::string line, type, from, from_orient, to;
    int count = 0;
    while (std::getline(in, line)) {
        std::istringstream fields(line);
        fields >> type >> from >> from_orient >> to;
        if (type != "L") continue;
        if (to_sink && to == "2") count += 1;
        if (!to_sink && from == "0") count += 1;
    }
    return count;
}

static theseus::MSAOptions make_options(theseus::EndCondition end_condition, bool reverse = false, bool density_drop = false) {
    theseus::MSAOptions options;
    options.end_condition = end_condition;
    options.reverse       = reverse;
    options.density_drop  = density_drop;
    return options;
}

TEST_CASE("Check MSA end conditions") {
    // The graph has 13 characters and the sequence lacks the last one
    const std::string initial_seq = "ACCCGTAAAAGGG";
    const std::string shorter_seq = "ACCCGTAAAAGG";
    theseus::Penalties penalties(0, 2, 3, 1);
    theseus::Heuristics heuristics;

    SUBCASE("Full query (default): the end of the graph is free") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        theseus::Alignment alignment = aligner.align(shorter_seq, 1, theseus::MSAOptions());
        CHECK(alignment.theseus_status == THESEUS_STATUS_ALG_COMPLETED);
        CHECK(alignment.score == 0);
    }

    SUBCASE("The default options are those of the call without options") {
        theseus::TheseusMSA aligner_1(penalties, heuristics, initial_seq, 1);
        theseus::TheseusMSA aligner_2(penalties, heuristics, initial_seq, 1);
        theseus::Alignment alignment_1 = aligner_1.align(shorter_seq, 1, false);
        theseus::Alignment alignment_2 = aligner_2.align(shorter_seq, 1, theseus::MSAOptions());
        CHECK(alignment_1.score == alignment_2.score);
        CHECK(alignment_1.edit_op == alignment_2.edit_op);
        CHECK(alignment_1.path == alignment_2.path);
    }

    SUBCASE("End to end: the end of the graph must be reached") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        theseus::Alignment alignment = aligner.align(shorter_seq, 1, make_options(theseus::EndCondition::EndToEnd));
        CHECK(alignment.theseus_status == THESEUS_STATUS_ALG_COMPLETED);
        CHECK(alignment.score == 4);                                     // Skip the last character: 3 + 1
        CHECK(alignment.path.back() == 2);                               // The path reaches the sink
    }

    SUBCASE("Ends free: the query is consumed or the end of the graph is reached") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        // The query is consumed before the end of the graph
        theseus::Alignment alignment = aligner.align(shorter_seq, 1, make_options(theseus::EndCondition::EndsFree));
        CHECK(alignment.score == 0);
        // The end of the graph is reached before the query is consumed
        alignment = aligner.align(initial_seq + "TTT", 1, make_options(theseus::EndCondition::EndsFree));
        CHECK(alignment.theseus_status == THESEUS_STATUS_ALG_COMPLETED);
        CHECK(alignment.score == 0);
        CHECK(alignment.path.back() == 2);
        // With the other end conditions the extra characters have to be aligned
        alignment = aligner.align(initial_seq + "TTT", 1, make_options(theseus::EndCondition::EndToEnd));
        CHECK(alignment.score == 6);                                     // 3 + 3
    }

    SUBCASE("Links to the sink of the POA graph") {
        const std::string prefix = "ACCCGT";    // First part of the graph
        // Full query: the sequence is linked to the sink from its last character
        theseus::TheseusMSA full_query(penalties, heuristics, initial_seq, 1);
        full_query.align(prefix, 1, make_options(theseus::EndCondition::FullQuery));
        CHECK(count_links(full_query, true) == 2);
        // End to end: the same, as it skips the rest of the graph
        theseus::TheseusMSA end_to_end(penalties, heuristics, initial_seq, 1);
        end_to_end.align(prefix, 1, make_options(theseus::EndCondition::EndToEnd));
        CHECK(count_links(end_to_end, true) == 2);
        // Ends free: there is no link, the sequence does not reach the end of the graph
        theseus::TheseusMSA ends_free(penalties, heuristics, initial_seq, 1);
        ends_free.align(prefix, 1, make_options(theseus::EndCondition::EndsFree));
        CHECK(count_links(ends_free, true) == 1);
        // In all the cases the MSA reproduces the sequence
        CHECK(ungapped_rows(full_query)[1] == prefix);
        CHECK(ungapped_rows(end_to_end)[1] == prefix);
        CHECK(ungapped_rows(ends_free)[1] == prefix);
    }
}

TEST_CASE("Check MSA reverse alignments") {
    // The sequences cover the last part of the graph
    const std::string initial_seq = "ACCCGTAAAAGGG";
    const std::string suffix      = "GTAAAAGGG";
    theseus::Penalties penalties(0, 2, 3, 1);
    theseus::Heuristics heuristics;

    SUBCASE("Full query: the start of the graph is free") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        theseus::Alignment alignment = aligner.align(suffix, 1, make_options(theseus::EndCondition::FullQuery, true));
        CHECK(alignment.theseus_status == THESEUS_STATUS_ALG_COMPLETED);
        CHECK(alignment.score == 0);
        CHECK(alignment.edit_op == std::vector<char>(9, 'M'));
        CHECK(alignment.path.back() == 2);                               // In the order of the graph
        CHECK(count_links(aligner, false) == 2);                         // Linked to the source
        CHECK(ungapped_rows(aligner)[1] == suffix);
    }

    SUBCASE("Ends free: the sequence is not linked to the source") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        theseus::Alignment alignment = aligner.align(suffix, 1, make_options(theseus::EndCondition::EndsFree, true));
        CHECK(alignment.score == 0);
        CHECK(count_links(aligner, false) == 1);
        CHECK(ungapped_rows(aligner)[1] == suffix);
    }

    SUBCASE("End to end: the start of the graph must be reached") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        theseus::Alignment alignment = aligner.align(suffix, 1, make_options(theseus::EndCondition::EndToEnd, true));
        CHECK(alignment.score == 7);                                     // Skip the first 4 characters: 3 + 4
        CHECK(alignment.path.front() == 0);                              // The path reaches the source
        CHECK(ungapped_rows(aligner)[1] == suffix);
    }

    SUBCASE("Forward and reverse sequences in the same MSA") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        aligner.align("ACCCGTAAA", 1, make_options(theseus::EndCondition::EndsFree));         // First part
        aligner.align(suffix, 1, make_options(theseus::EndCondition::EndsFree, true));        // Last part
        aligner.align("ACCCGTCAAAGGG", 1, make_options(theseus::EndCondition::EndToEnd));     // Full, with a mismatch
        std::vector<std::string> rows = ungapped_rows(aligner);
        REQUIRE(rows.size() == 5);                                                            // + consensus
        CHECK(rows[1] == "ACCCGTAAA");
        CHECK(rows[2] == suffix);
        CHECK(rows[3] == "ACCCGTCAAAGGG");
    }
}

TEST_CASE("Check MSA density drop") {
    // Pseudo-random bases (fixed seed)
    unsigned int state = 12345;
    std::string initial_seq, junk;
    for (int i = 0; i < 250; ++i) { state = state * 1103515245u + 12345u; initial_seq += "ACGT"[(state >> 16) & 3]; }
    for (int i = 0; i < 200; ++i) { state = state * 1103515245u + 12345u; junk += "ACGT"[(state >> 16) & 3]; }
    theseus::Penalties penalties(0, 2, 3, 1);
    theseus::Heuristics heuristics;

    SUBCASE("Without the heuristic the whole sequence is aligned") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        theseus::Alignment alignment = aligner.align(initial_seq + junk, 1, make_options(theseus::EndCondition::EndToEnd));
        CHECK(alignment.theseus_status == THESEUS_STATUS_ALG_COMPLETED);
        CHECK(alignment.score == 3 + 200);                               // One insertion of the junk
    }

    SUBCASE("A dropped alignment is partial, and it is added to the POA graph") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        const std::string seq = initial_seq + junk;
        theseus::Alignment alignment = aligner.align(seq, 1, make_options(theseus::EndCondition::EndToEnd, false, true));
        CHECK(alignment.theseus_status == THESEUS_STATUS_END_UNREACHABLE);
        // The score is the one of the partial alignment (not of the whole sequence)
        CHECK(alignment.score < 3 + 200);
        // Only the aligned part of the sequence is in the MSA
        std::vector<std::string> rows = ungapped_rows(aligner);
        REQUIRE(rows.size() == 3);                                       // + consensus
        CHECK(rows[1].size() > 0);
        CHECK(rows[1].size() < seq.size());
        CHECK(seq.compare(0, rows[1].size(), rows[1]) == 0);
        // The following sequences can be added with the POA graph in sync
        aligner.align(initial_seq, 1);
        rows = ungapped_rows(aligner);
        REQUIRE(rows.size() == 4);
        CHECK(rows[2] == initial_seq);
    }

    SUBCASE("Reverse alignment: the first part of the sequence is the junk") {
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);
        const std::string seq = junk + initial_seq;
        theseus::Alignment alignment = aligner.align(seq, 1, make_options(theseus::EndCondition::EndToEnd, true, true));
        CHECK(alignment.theseus_status == THESEUS_STATUS_END_UNREACHABLE);
        std::vector<std::string> rows = ungapped_rows(aligner);
        REQUIRE(rows.size() == 3);
        CHECK(rows[1].size() > 0);
        CHECK(rows[1].size() < seq.size());
        // The aligned part is the end of the sequence
        CHECK(seq.compare(seq.size() - rows[1].size(), rows[1].size(), rows[1]) == 0);
    }
}
