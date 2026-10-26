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
        theseus::Alignment alignment = aligner.align(new_seq, false);

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
        theseus::Alignment alignment = aligner.align(seq_1, false);
        CHECK(alignment.score == 2); // Check score

        // Align and check sequence 2
        alignment = aligner.align(seq_2, false);
        CHECK(alignment.score == 6); // Check score

        // Align and check sequence 3
        alignment = aligner.align(seq_3, false);
        CHECK(alignment.score == 0); // Check score

        // Align and check sequence 4
        alignment = aligner.align(seq_4, false);
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
        // alignment = aligner.align(seq_4, false);
        // CHECK(alignment.score == 9); // Check score
    }
}