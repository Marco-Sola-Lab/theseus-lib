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
 * @brief Parse the header of a sequence in the Pericles format:
 *          "> reverse end_condition weight density_drop lag_pruning"
 * The five integers must be present. The flags (reverse, density_drop and
 * lag_pruning) must be 0 or 1, the end condition 0 (full query), 1 (ends free) or
 * 2 (end to end), and the weight positive. Otherwise it is not a Pericles header
 * (e.g. a sequence name) and the function returns false.
 *
 * @param line     Header line of the fasta file (starts with '>')
 * @param weight   Weight of the sequence
 * @param options  Options of the alignment (reverse, end condition, density drop and lag pruning)
 */
bool parse_header(
    const std::string &line,
    int &weight,
    theseus::MSAOptions &options)
{
    std::istringstream iss(line);
    char header;
    int reverse, end_condition, density_drop, lag_pruning;
    if (!(iss >> header >> reverse >> end_condition >> weight >> density_drop >> lag_pruning)) return false;
    const bool valid_flags = (reverse == 0 || reverse == 1) && (end_condition >= 0 && end_condition <= 2) &&
                             (density_drop == 0 || density_drop == 1) && (lag_pruning == 0 || lag_pruning == 1);
    if (!valid_flags || weight < 1) return false;
    options.reverse       = (reverse == 1);
    // End condition of the header: 0 full query, 1 ends free, 2 end to end
    if (end_condition == 0) options.end_condition = theseus::EndCondition::FullQuery;
    else if (end_condition == 1) options.end_condition = theseus::EndCondition::EndsFree;
    else options.end_condition = theseus::EndCondition::EndToEnd;
    options.density_drop  = (density_drop == 1);
    options.lag_pruning   = (lag_pruning == 1);
    return true;
}


/**
 * @brief Read the sequences from a file. If the header of a sequence is in the
 * Pericles format (see parse_header) its weight and options are read from it.
 * Otherwise the sequence has weight 1 and the default options (the lag pruning
 * is the one of the command line).
 *
 * @param sequences Vector to store the sequences
 * @param weights Vector to store the weights of the sequences
 * @param options Vector to store the alignment options of the sequences
 * @param args Arguments containing the file name
 */
void read_sequences(
    std::vector<std::string> &sequences,
    std::vector<int> &weights,
    std::vector<theseus::MSAOptions> &options,
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
            // Default values, replaced by the ones of the header (if it has them)
            int weight = 1;
            theseus::MSAOptions seq_options;
            seq_options.lag_pruning = args.lag_pruning;
            if (!parse_header(line, weight, seq_options)) {
                weight = 1;
                seq_options = theseus::MSAOptions();
                seq_options.lag_pruning = args.lag_pruning;
            }
            weights.push_back(weight);
            options.push_back(seq_options);
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
    std::cout << "Usage: theseus_msa [OPTIONS]\n"
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

                 " Per-sequence options (optional, in the header of each sequence of the fasta file):\n"
                 "  > <reverse> <end_condition> <weight> <density_drop> <lag_pruning>\n"
                 "  If the header has these five integers they set the options of the sequence: reverse\n"
                 "  alignment (0 or 1), end condition (0: full query, 1: ends free, 2: end to end), weight,\n"
                 "  density drop (0 or 1) and lag pruning (0 or 1). Otherwise: weight 1, forward alignment\n"
                 "  that ends when the whole sequence is consumed (full query), no density drop, and the\n"
                 "  lag pruning of the command line.\n"

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
    std::vector<theseus::MSAOptions> options;
    read_sequences(sequences, weights, options, args);
    if (sequences.empty()) {
        std::cerr << "No sequences found in the dataset file\n";
        return 1;
    }

    // Prepare the data
    std::vector<theseus::Alignment> alignments(sequences.size());
    std::string_view initial_seq = sequences[0];
    // The options of the first sequence are not used: it only creates the POA graph
    theseus::TheseusMSA aligner(penalties, heuristics, initial_seq, weights[0]);

    // Alignment with Theseus
    for (int j = 1; j < sequences.size(); ++j) {
        std::cout << "Processing sequence " << j << std::endl;
        alignments[j] = aligner.align(sequences[j], weights[j], options[j]);
        if (alignments[j].theseus_status == THESEUS_STATUS_END_UNREACHABLE) {
            std::cerr << "Alignment " << j << " dropped by the density heuristic (partial alignment added)" << std::endl;
        }
        else if (alignments[j].theseus_status != THESEUS_STATUS_ALG_COMPLETED) {
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