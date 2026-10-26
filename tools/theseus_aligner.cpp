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
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <cassert>
#include <unordered_map>


#include "theseus/alignment.h"
#include "theseus/graph.h"
#include "theseus/penalties.h"
#include "theseus/heuristics.h"
#include "theseus/theseus_aligner.h"


// Control of the output.
#define AVOID_DP 0
#define AVOID_THESEUS 0
#define PRINT_ALIGNMENTS 0

using NodeId = theseus::Graph::NodeId;

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
    bool density_drop = false;
    bool lag_pruning  = false;
    // I/O
    std::string graph_file;
    std::string sequences_and_positions_file;
    std::string output_file;
};


NodeId node_name_to_id(theseus::Graph &graph,
                       const std::string &name,
                       const std::string &dna_seq,
                       std::unordered_map<std::string, NodeId> &name_to_id,
                       std::unordered_map<NodeId, std::string> &node_names)
{
    auto seq_ptr = name_to_id.find(name);
    // If the name is not found, add it
    if (seq_ptr == name_to_id.end())
    {
        // Check that lengths are consistent
        assert(name_to_id.size() == graph.nnodes());
        // Add the new node to the graph
        NodeId id = graph.add_node(dna_seq);
        // Add the new name
        name_to_id[name] = id;
        node_names[id]   = name;
        return id;
    }

    // Check consistency
    assert(seq_ptr->second < graph.nnodes());
    assert(name == node_names[seq_ptr->second]);
    return seq_ptr->second; // Second is the value of the key-value pair
}


/**
 * @brief Return the reverse complement of a DNA string.
 *
 * @param dna_string The DNA string to reverse complement.
 * @return std::string The reverse complement of the input DNA string.
 */
std::string reverse_complement(const std::string &dna_string) {
    std::string rev_comp = dna_string;
    std::reverse(rev_comp.begin(), rev_comp.end()); 	// Reverse the DNA string
    for (char &c : rev_comp) {							// Complement the DNA string
        switch (c) {
            case 'A': case 'a': c = 'T'; break;
            case 'T': case 't': c = 'A'; break;
            case 'C': case 'c': c = 'G'; break;
            case 'G': case 'g': c = 'C'; break;
        }
    }
    return rev_comp;
}


// Construct a graph object from a  GFA file stream
void graph_from_gfa_stream(std::istream &gfa_stream,
                                     theseus::Graph &graph,
                                     std::unordered_map<std::string, NodeId> &name_to_id,
                                     std::unordered_map<NodeId, std::string> &node_names) {
    std::string line;
    while (gfa_stream.good()) {
        std::getline(gfa_stream, line);
        if (line.size() == 0 && !gfa_stream.good())
            break;

        // Only Segments and Links are supported
        if (line.size() == 0 || (line[0] != 'S' && line[0] != 'L'))
            continue;

        // Parse segment data
        if (line[0] == 'S')
        {
            std::stringstream sstr{line};
            std::string type, name, dna_seq;
            sstr >> type;
            assert(type == "S");
            sstr >> name >> dna_seq;
            // Consider first the forward orientation
            name = name + "+";
            NodeId id = node_name_to_id(graph, name, dna_seq, name_to_id, node_names);
            // Warning on empty nodes
            if (dna_seq == "*")
                std::cerr << std::string{"Nodes without sequence (*) are not currently supported (nodeid " + std::to_string(id) + ")"};
            assert(dna_seq.size() >= 1);
            // We add the reverse orientation
            std::string rev_name = name.substr(0, name.size() - 1) + "-";
            std::string rev_dna_seq = reverse_complement(dna_seq);
            NodeId rev_id = node_name_to_id(graph, rev_name, rev_dna_seq, name_to_id, node_names);
        }
        // Parse link data. We add both edges (fromstr+fromstart, tostr+toend)
        // and (tostr+toend, fromstr+fromstart), to support bidirectedness
        if (line[0] == 'L')
        {
            std::stringstream sstr{line};
            std::string type, fromstr_forward, tostr_forward, fromstr_reverse, tostr_reverse, fromstart, toend, overlapstr;
            int overlap = 0;
            sstr >> type;
            sstr >> fromstr_forward >> fromstart >> tostr_forward >> toend >> overlapstr;
            // Assess if the read data is consistent with the format
            assert(type == "L");
            assert(fromstart == "+" || fromstart == "-");
            assert(toend == "+" || toend == "-");
            // Set name ids
            fromstr_forward = fromstr_forward + fromstart;
            tostr_forward   = tostr_forward + toend;
            fromstr_reverse = fromstr_forward.substr(0, fromstr_forward.size() - 1) + (fromstart == "+" ? "-" : "+");
            tostr_reverse   = tostr_forward.substr(0, tostr_forward.size() - 1) + (toend == "+" ? "-" : "+");
            // Get the node ids for the forward and reverse orientations of the edge
            NodeId from_forward = node_name_to_id(graph, fromstr_forward, "", name_to_id, node_names);
            NodeId to_forward   = node_name_to_id(graph, tostr_forward  , "", name_to_id, node_names);
            NodeId from_reverse = node_name_to_id(graph, fromstr_reverse, "", name_to_id, node_names);
            NodeId to_reverse   = node_name_to_id(graph, tostr_reverse  , "", name_to_id, node_names);
            // Check overlap (currently only accept 0M)
            if (overlapstr != "0M") {
                std::cerr << "Currently, only edge overlaps of 0M are supported (non supported overlap: " + overlapstr + ")" << std::endl;
            }
            // Store the edges
            // Forward edge
            int frompos = (int)from_forward;
            int topos   = (int)to_forward;
            graph.add_edge(from_forward, to_forward);
            // Reverse edge
            frompos = (int)to_reverse;
            topos   = (int)from_reverse;
            graph.add_edge(to_reverse, from_reverse);
        }
    }
    // Check that nodes are not empty
    for (NodeId id : graph.nodes())
    {
        auto node = graph.node(id);
        if (node.sequence.size() > 0)
            continue;
        std::string name = node_names[id];
        if (name.back() == '+') {
            std::cerr << std::string{"Node " + name + " is present in edges but missing in nodes"} << std::endl;
        }
    }
    // Validate that all edges connect existing nodes TODO:
}


// Read sequence data
void read_seq_pos_data(
    std::ifstream &sp_file,
    std::vector<std::string> &sequences,
    std::vector<NodeId> &start_nodes,
    std::vector<int> &start_offsets,
    std::unordered_map<std::string, NodeId> &name_to_id)
{
    if (!sp_file.is_open()) {
        std::cerr << "Could not open dataset file\n";
        return;
    }

    std::string sequence, line; // Value and metadata of the sequence
    std::istringstream iss(line);
    std::string vertex, orientation;
    int offset;

    // Read sequences
    int num = 0;
    while (getline(sp_file, line))
    {
        if (line.empty())
            continue;

        if (line[0] == '>')
        {
            // Read positional data
            iss.clear();
            iss.str(line.substr(1)); // Skip the '>'
            if (!(iss >> vertex >> offset >> orientation)) {
                std::cerr << "Error reading position line: " << line << std::endl;
                continue;
            }
            if (orientation != "+" && orientation != "-") {
                std::cerr << "Invalid orientation in line: " << line << std::endl;
                continue;
            }

            vertex = vertex + orientation; // We store the orientation in the vertex name);
            NodeId node_id = name_to_id[vertex];
            start_nodes.push_back(node_id);
            start_offsets.push_back(offset);

            // Store previous sequence if any
            if (num > 0) sequences.push_back(sequence);
            sequence.clear();
            num += 1;
        }
        else
        {
            sequence += line;   // The sequnce may span several lines
        }
    }

    // Store the last sequence
    if (num > 0) {
      sequences.push_back(sequence);
    }

    // Close the file
    sp_file.close();
}


/**
 * @brief Print the help message.
 */
void help() {
    std::cout << "Usage: theseus_aligner [OPTIONS]\n"
                 "Options:\n"
                 "  Penalties:\n"
                 "  -m, --match <int>            The match penalty                                       [default=0]\n"
                 "  -x, --mismatch <int>         The mismatch penalty                                    [default=2]\n"
                 "  -o, --gapo <int>             The gap open penalty                                    [default=3]\n"
                 "  -e, --gape <int>             The gap extension penalty                               [default=1]\n"
                 "  -O, --gapo2 <int>            The second gap open penalty (dual affine)               [default=none]\n"
                 "  -E, --gape2 <int>            The second gap extension penalty (dual affine)          [default=none]\n\n"

                 "  I/O:\n"
                 "  -g, --graph_file <file>      Graph file in .gfa format                               [Required]\n"
                 "  -s, --sequences_file <file>  Sequences and starting positons in .fasta format        [Required]\n"
                 "  -f, --output_file <file>     Output file                                             [Required]\n\n"

                 " Heuristics:\n"
                 "  -d  --density_heuristic     Activate the drop heuristic based on advancement density.            \n"
                 "  -l  --lag_pruning           Activate the pruning of diagonals lagging behind int the alignment.  \n";
}

CMDArgs parse_args(int argc, char *const *argv) {
    static const option long_options[] = {{"match", required_argument, 0, 'm'},
                                          {"mismatch", required_argument, 0, 'x'},
                                          {"gapo", required_argument, 0, 'o'},
                                          {"gape", required_argument, 0, 'e'},
                                          {"gapo2", required_argument, 0, 'O'},
                                          {"gape2", required_argument, 0, 'E'},
                                          {"graph_file", required_argument, 0, 'g'},
                                          {"sequences_file", required_argument, 0, 's'},
                                          {"output_file", required_argument, 0, 'f'},
                                          {"lag_pruning", no_argument, 0, 'l'},
                                          {"density_heuristic", no_argument, 0, 'd'},
                                          {0, 0, 0, 0}};

    CMDArgs args;

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "m:x:o:e:O:E:g:s:f:ld", long_options, &option_index)) != -1) {
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
            case 'g':
                args.graph_file = optarg;
                break;
            case 's':
                args.sequences_and_positions_file = optarg;
                break;
            case 'f':
                args.output_file = optarg;
                break;
            case 'l':
                args.lag_pruning = true;
                break;
            case 'd':
                args.density_drop = true;
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

    if (args.graph_file.empty() || args.sequences_and_positions_file.empty() || args.output_file.empty()) {
        std::cerr << "Missing required arguments\n";
        help();
        return 1;
    }

    // Parse penalties
    // Both pieces of the second gap function are needed for dual affine gaps
    if ((args.gapo2 < 0) != (args.gape2 < 0)) {
        std::cerr << "The second gap open (-O) and extension (-E) penalties must be set together\n";
        return 1;
    }
    const bool dual_affine = args.gapo2 >= 0;
    theseus::Penalties penalties = dual_affine
        ? theseus::Penalties(args.match, args.mismatch, args.gapo, args.gape, args.gapo2, args.gape2)
        : theseus::Penalties(args.match, args.mismatch, args.gapo, args.gape);
    // Parse heuristics
    theseus::Heuristics heuristics;
    // Manage input/output files
    std::ifstream graph_file(args.graph_file);
    std::ifstream sp_file(args.sequences_and_positions_file);
    std::ofstream output_file(args.output_file);
    // Construct the graph
    theseus::Graph graph;
    std::unordered_map<std::string, NodeId> name_to_id;
    std::unordered_map<NodeId, std::string> node_names;
    graph_from_gfa_stream(graph_file, graph, name_to_id, node_names);
    // Prepare the aligner
    theseus::TheseusAligner aligner(penalties, heuristics, std::move(graph));
    // Read queries data
    std::vector<std::string> sequences;
    std::vector<NodeId> start_nodes;
    std::vector<int> start_offsets;
    read_seq_pos_data(sp_file, sequences, start_nodes, start_offsets, name_to_id);
    // Align the sequences
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    int num_sequences = sequences.size();
    theseus::Alignment alignment;
    for (int i = 0; i < num_sequences; ++i) {
        // Perform alignment
        std::cout << "Seq " << i << std::endl;
        alignment = aligner.align(sequences[i], start_nodes[i], start_offsets[i], args.density_drop, args.lag_pruning);
        std::cout << "Score = " << alignment.score << std::endl << std::endl;
        if (alignment.theseus_status != THESEUS_STATUS_ALG_COMPLETED) {
            std::cerr << "Alignment " << i << " with status " << alignment.theseus_status << " did not complete successfully" << std::endl;
        }
        aligner.print_alignment_as_gaf(alignment, output_file, "seq_" + std::to_string(i), node_names);
    }
    // End time measurement
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;
    // Close files
    graph_file.close();
    sp_file.close();
    output_file.close();
    return 0;
}