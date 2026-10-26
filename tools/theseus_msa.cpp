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


#include <getopt.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <sstream>

#include "theseus/alignment.h"
#include "theseus/heuristics.h"
#include "theseus/penalties.h"
#include "theseus/theseus_msa_aligner.h"

#include <vector>

// Command line arguments
struct CMDArgs {
    // Penalties
    int match = 0;
    int mismatch = 2;
    int gapo = 3;
    int gape = 1;
    // Second piece of the gap function (dual affine). Both must be set to use it
    int gapo2 = -1;
    int gape2 = -1;
    // Heuristics
    bool lag_pruning = false;
    // I/O
    int output_type = 0;        // 0: MSA, 1: GFA, 2: Consensus, 3: Dot
    std::string sequences_file;
    std::string output_file;
};


/**
 * @brief Read the sequences from a file.
 *
 * @param sequences Vector to store the sequences
 * @param args Arguments containing the file name
 */
void read_sequences(
    std::vector<std::string> &sequences,
    std::vector<int> &weights,
    CMDArgs &args)
{

    // Open the file containing the sequences
    std::ifstream sequences_file(args.sequences_file);

    if (!sequences_file.is_open()) {
        std::cerr << "Could not open dataset file\n";
        return;
    }

    std::string sequence, line; // Value and metadata of the sequence

    // TODO: Allow for several sequence formats
    int num = 0;
    while (getline(sequences_file, line))
    {
        if (line.empty())
            continue;

        if (line[0] == '>')
        {
            if (num > 0) sequences.push_back(sequence);
            sequence.clear();
            num += 1;
            // std::istringstream iss(line);
            // char header;
            // int weight;
            // iss >> header >> weight;
            // weights.push_back(weight);
        }
        else
        {
            sequence += line;   // The sequnce may span several lines
        }
    }

    // Store the last sequencealignments
    if (num > 0) {
      sequences.push_back(sequence);
    }

    // Close the file
    sequences_file.close();
}


/**
 * @brief Print the help message.
 */
void help() {
    std::cout << "Usage: pericles [OPTIONS]\n"
                 "Options:\n"
                 "  -m, --match <int>           The match penalty                                       [default=0]\n"
                 "  -x, --mismatch <int>        The mismatch penalty                                    [default=2]\n"
                 "  -o, --gapo <int>            The gap open penalty                                    [default=3]\n"
                 "  -e, --gape <int>            The gap extension penalty                               [default=1]\n"
                 "  -O, --gapo2 <int>           The second gap open penalty (dual affine)               [default=none]\n"
                 "  -E, --gape2 <int>           The second gap extension penalty (dual affine)          [default=none]\n"
                 "  -t, --output_type <int>     The output format of the multiple alignment             [default=0=MSA]\n"
                 "                               0: MSA: Standard Multiple Sequence Alignment format,\n"
                 "                               1: GFA: Output the resulting POA graph in GFA format,\n"
                 "                               2: Consensus - Heaviest Bundle: Output the consensus sequence \n"
                 "                                  based on the heaviest bundle algorithm,\n"
                 "                               3: Consensus - Weighted Majority Voting: Output the consensus \n"
                 "                                  sequence based on the weighted majority voting algorithm. \n"
                 "                                  Moreover, it prints each column's weight and the gapped \n"
                 "                                  consensus,\n"
                 "                               4: Dot: Output in .dot format for visualization purposes.\n"
                 "                                       Only tractable for small graphs\n"
                 "  -f, --output <file>         Output file                                             [Required]\n"
                 "  -s, --sequences <file>      Dataset file                                            [Required]\n\n"

                 " Heuristics:\n"
                 "  -l  --lag_pruning           Activate the pruning of diagonals lagging behind int the alignment.  \n";
}

CMDArgs parse_args(int argc, char *const *argv) {
    static const option long_options[] = {{"match", required_argument, 0, 'm'},
                                          {"mismatch", required_argument, 0, 'x'},
                                          {"gapo", required_argument, 0, 'o'},
                                          {"gape", required_argument, 0, 'e'},
                                          {"gapo2", required_argument, 0, 'O'},
                                          {"gape2", required_argument, 0, 'E'},
                                          {"output_type", required_argument, 0, 't'},
                                          {"sequences", required_argument, 0, 's'},
                                          {"output", required_argument, 0, 'f'},
                                          {"lag_pruning", no_argument, 0, 'l'},
                                          {0, 0, 0, 0}};

    CMDArgs args;

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "m:x:o:e:O:E:t:s:f:l", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'm':
                args.match = std::stoi(optarg);
                break;
            case 'x':
                args.mismatch = std::stoi(optarg);
                break;
            case 'o':
                args.gapo = std::stoi(optarg);
                break;
            case 'e':
                args.gape = std::stoi(optarg);
                break;
            case 'O':
                args.gapo2 = std::stoi(optarg);
                break;
            case 'E':
                args.gape2 = std::stoi(optarg);
                break;
            case 't':
                args.output_type = std::stoi(optarg);
                break;
            case 's':
                args.sequences_file = optarg;
                break;
            case 'f':
                args.output_file = optarg;
                break;
            case 'l':
                args.lag_pruning = true;
                break;
            default:
                std::cerr << "Invalid option" << std::endl;
                exit(1);
        }
    }

    return args;
}


int main(int argc, char *const *argv) {
    // Parsing
    CMDArgs args = parse_args(argc, argv);

    if (args.sequences_file.empty() || args.output_file.empty()) {
        std::cerr << "Missing required arguments\n";
        help();
        return 1;
    }

    // Check ouput type is correct
    if (args.output_type < 0 || args.output_type > 4) {
        std::cerr << "Output type must be 0 (MSA), 1 (GFA), 2 (Consensus), 3 (Majority Voting Consensus) or 4 (Dot)\n";
        return 1;
    }

    // Define alignment penalties
    // Both pieces of the second gap function are needed for dual affine gaps
    if ((args.gapo2 < 0) != (args.gape2 < 0)) {
        std::cerr << "The second gap open (-O) and extension (-E) penalties must be set together\n";
        return 1;
    }
    const bool dual_affine = args.gapo2 >= 0;
    theseus::Penalties penalties = dual_affine
        ? theseus::Penalties(args.match, args.mismatch, args.gapo, args.gape, args.gapo2, args.gape2)
        : theseus::Penalties(args.match, args.mismatch, args.gapo, args.gape);
    // Determine heuristics
    theseus::Heuristics heuristics;
    // Read the sequences for the MSA
    std::vector<std::string> sequences;
    std::vector<int> weights;
    read_sequences(sequences, weights, args);

    // Prepare the data
    std::vector<theseus::Alignment> alignments(sequences.size());
    std::string_view initial_seq = sequences[0];
    theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, 1);

    // Alignment with Theseus
    for (int j = 1; j < sequences.size(); ++j) {
        std::cout << "Processing sequence " << j << std::endl;
        alignments[j] = aligner.align(sequences[j], 1, args.lag_pruning);
        if (alignments[j].theseus_status != THESEUS_STATUS_ALG_COMPLETED) {
            std::cerr << "Alignment " << j << " with status " << alignments[j].theseus_status << " did not complete successfully" << std::endl;
        }
        std::cout << "Score = " << alignments[j].score << std::endl << std::endl;
    }

    // Print the output
    std::ofstream output_file(args.output_file);
    if (args.output_type == 0) {
        aligner.print_as_msa(output_file);
    } else if (args.output_type == 1) {
        aligner.print_as_gfa(output_file);
    } else if (args.output_type == 2) {
        std::string consensus = aligner.heaviest_bundle_consensus();
        output_file << ">Consensus\n" << consensus << "\n";
    } else if (args.output_type == 3) {
        std::vector<int> consensus_weights;
        std::string consensus_sequence, consensus_sequence_gapped;
        aligner.majority_voting_consensus(consensus_weights, consensus_sequence, consensus_sequence_gapped);
        output_file << ">Consensus\n" << consensus_sequence << "\n";
        output_file << ">Gapped_consensus\n" << consensus_sequence_gapped << "\n";
        output_file << ">Consensus_weights\n";
        for (size_t i = 0; i < consensus_weights.size(); ++i) {
            output_file << consensus_weights[i];
            if (i < consensus_weights.size() - 1) {
                output_file << " ";
            }
        }
        output_file << "\n";

    } else if (args.output_type == 4) {
        aligner.print_as_dot(output_file);
    }

    return 0;
}