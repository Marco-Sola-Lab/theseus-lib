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


#pragma once

#include <memory>
#include <string>
#include <queue>
#include <set>
#include <algorithm>
#include <ranges>
#include <stack>
#include <utility>
#include <tuple>

#include "theseus/alignment.h"
#include "theseus/penalties.h"
#include "theseus/heuristics.h"
#include "theseus/graph.h"

#include "beyond_scope.h"
#include "cell.h"
#include "scope.h"
#include "scratchpad.h"
#include "vertices_data.h"
#include "wavefront.h"
#include "internal_penalties.h"
#include "msa.h"

namespace theseus {

using NodeId       = Graph::NodeId;
using NodeView     = Graph::NodeView;
using NodeIdRange  = Graph::NodeIdRange;
using SequenceView = Graph::SequenceView;

class TheseusAlignerImpl {
public:
    TheseusAlignerImpl(const Penalties &penalties,
                       const Heuristics &heuristics,
                       Graph &&graph,
                       int  initial_weight,
                       bool is_msa);

    /**
     * @brief Main alignment function. Aligns the given sequence to the graph
     * starting at the specified node and offset.
     *
     * @param seq                Sequence to be aligned
     * @param weight             Weight of the sequence to be aligned (used for MSA)
     * @param reverse_alignment  Whether to perform reverse alignment
     * @param is_ends_free       Whether to allow a free end on the "end" of the graph
     *
     * @return                  Alignment object
     */
    Alignment align(std::string_view seq,
                    // Seq-to-graph parameters
                    int  start_node,
                    int  start_offset,
                    // MSA parameters
                    int  weight = 1,
                    bool is_ends_free = false,
                    // Common parameters
                    bool reverse_alignment = false,
                    bool density_drop_active = false,
                    bool lag_pruning_active = false);

    /**
     * @brief Output the current graph in GFA format.
     *
     * @param gfa_output  Output stream to write the graph in GFA format
     * @param node_names  Vector containing the names of the nodes in the graph (in the same order as their IDs)
     */
    void print_as_gfa(std::ostream &gfa_output);

    /**
     * @brief Output the current MSA in MSA format.
     *
     * @param out_stream  Output stream to write the MSA in MSA format
     */
    void print_as_msa(std::ostream &out_stream);

    /**
     * @brief Return the consensus sequence from the current MSA.
     *
     * @return std::string Consensus sequence
     */
    std::string heaviest_bundle_consensus();

    /**
     * @brief Get the weighted majority voting consensus sequence from the current
     * MSA.
     *
     */
    void majority_voting_consensus(std::vector<int> &consensus_weights,
                                   std::string &consensus_sequence,
                                   std::string &consensus_sequence_gapped);

    /**
     * @brief Print the graph in dot (graphviz) format
     *
     */
    void print_code_graphviz(std::ostream &out_stream);

    /**
     * @brief Print the resulting alignment in GAF format.
     *
     * @param alignment Alignment to be printed
     * @param out_stream Output stream where the alignment will be printed
     */
    void print_as_gaf(
            theseus::Alignment &alignment,
            std::ostream &out_stream,
            std::string seq_name,
            std::unordered_map<NodeId, std::string> &node_names);

private:
    // Pending I jump to a neighbour (used to cross empty vertices without recursion)
    struct PendingIJump {
        NodeId       node_id;
        Cell::pos_t  prev_pos;
        Cell::Matrix from_matrix;
    };

    // Pending diagonal to extend (stack of extend_diagonal)
    struct ExtendState {
        NodeId       node_id;
        Cell::pos_t  pos;
        Cell::Matrix from_matrix;
    };

    /**
     * @brief Initialize the data for a new alignment.
     *
     */
    void new_alignment(SequenceView seq,
                       bool reverse_alignment,
                       bool density_drop_active,
                       bool lag_pruning_active);

    /**
     * @brief Process a given vertex at a given _score. This means performing
     * the next and extend operations.
     *
     * @param curr_node_id
     */
    void process_vertex(
        NodeId curr_node_id);

    /**
     * @brief Compute the wave for a given score for all active vertices.
     *
     */
    void compute_new_wave();

    /**
     * @brief Sparsify the M data. This means storing the data in the scratchpad
     * to be later processed.
     *
     * @param dense_wf
     * @param offset_increase
     * @param shift_factor
     * @param cells_range
     * @param m
     * @param upper_bound
     */
    void sparsify_M_data(
        Cell::CellVector &dense_wf,
        int offset_increase,
        int shift_factor,
        Scope::range cells_range,
        int m,
        int upper_bound);

    /**
     * @brief Sparsify the jumps data. This means storing the data in the scratchpad
     * to be later processed.
     *
     * @param dense_wf
     * @param jumps_positions
     * @param offset_increase
     * @param shift_factor
     * @param m
     * @param upper_bound
     * @param from_matrix
     */
    void sparsify_jumps_data(
        Cell::CellVector &dense_wf,
        std::vector<Cell::pos_t> &jumps_positions,
        int offset_increase,
        int shift_factor,
        int m,
        int upper_bound,
        Cell::Matrix from_matrix);

    /**
     * @brief Sparsify the indel (coming from I or D) data. This means storing
     * the data in the scratchpad to be later processed.
     *
     * @param dense_wf
     * @param offset_increase
     * @param shift_factor
     * @param cells_range
     * @param m
     * @param upper_bound
     */
    void sparsify_indel_data(
        Cell::CellVector &dense_wf,
        int offset_increase,
        int shift_factor,
        Scope::range cells_range,
        int m,
        int upper_bound);

    /**
     * @brief Compute the next I matrix for a vertex v. This implies both sparsifying
     * the data in the scratchpad and storing it back on the new wavefront, once the
     * corresponding maximums and checks have been done.
     *
     * @param upper_bound // Maximum value of the diagonal
     * @param curr_node_id
     */
    template <Cell::Matrix matrix>
    void next_I(
        int upper_bound,
        NodeId curr_node_id);

    /**
     * @brief Compute the next D matrix for a vertex v. This implies both sparsifying
     * the data in the scratchpad and storing it back on the new wavefront, once the
     * corresponding maximums and checks have been done.
     *
     * @param upper_bound // Maximum value of the diagonal
     * @param curr_node_id
     */
    template <Cell::Matrix matrix>
    void next_D(
        int upper_bound,
        NodeId curr_node_id);

    // Compile-time dispatch of the data of the gap matrices (I, I2, D and D2)
    template <Cell::Matrix matrix>
    Cell::CellVector &gap_wf(int score);          // Scope wavefront
    template <Cell::Matrix matrix>
    Scope::RangeVector &gap_pos(int score);       // Scope positions
    template <Cell::Matrix matrix>
    Cell::CellVector &gap_jumps_wf();             // Jumps wavefront (I and I2)
    template <Cell::Matrix matrix>
    std::vector<std::vector<Cell::pos_t>> &gap_jumps_positions(NodeId node_id); // Jumps positions (I and I2)

    /**
     * @brief Compute the next M matrix for a vertex v. This implies both sparsifying
     * the data in the scratchpad and storing it back on the new wavefront, once the
     * corresponding maximums and checks have been done.
     *
     * @param upper_bound // Maximum value of the diagonal
     * @param curr_node_id
     */
    void next_M(
        int upper_bound,
        NodeId curr_node_id);

    /**
     * @brief Invalidate the diagonal associated to a jump in M, activate the newly
     * discovered vertices and store the jump in the neighbours.
     *
     * @param curr_node
     * @param jump_cell  Cell that reached the end of curr_node (passed by value)
     */
    void store_M_jump(
        NodeView curr_node,
        Cell jump_cell);

    /**
     * @brief Store an I (or I2) jump in a destination vertex. If the vertex is
     * empty, the jump continues: its neighbours are pushed to pending.
     *
     * @tparam matrix  Cell::Matrix::I or Cell::Matrix::I2
     * @param dest_id  Destination vertex
     * @param prev_pos  Position of the origin cell of the gap (for the backtrace)
     * @param from_matrix  Matrix of the origin cell of the gap
     * @param jump_cell  Jump cell, modified in place
     * @param pos_score  Position of the current score in the vertices data
     * @param pending  Pending jumps through empty vertices
     */
    template <Cell::Matrix matrix>
    void store_I_jump_in_vertex(
        NodeId dest_id,
        Cell::pos_t prev_pos,
        Cell::Matrix from_matrix,
        Cell &jump_cell,
        int pos_score,
        std::vector<PendingIJump> &pending);

    /**
     * @brief Invalidate the diagonal associated to a jump in I (or I2), activate
     * the newly discovered vertices and store the jump in the neighbours.
     *
     * @tparam matrix  Cell::Matrix::I or Cell::Matrix::I2
     * @param curr_node
     * @param jump_cell  Cell that reached the end of curr_node (passed by value)
     */
    template <Cell::Matrix matrix>
    void store_I_jump(
        NodeView curr_node,
        Cell jump_cell);

    /**
     * @brief Check and store I jumps (that is, those diagonals that have reached
     * the last column of a vertex for matrix I or I2).
     *
     * @tparam matrix  Cell::Matrix::I or Cell::Matrix::I2
     * @param curr_node
     * @param curr_wavefront
     * @param cell_range
     */
    template <Cell::Matrix matrix>
    void check_and_store_jumps(
        NodeView curr_node,
        Cell::CellVector &curr_wavefront,
        Scope::range cell_range);

    /**
     * @brief Longest Common Prefix of two sequences. The first sequence is always
     * the query _seq.
     *
     * @param text
     * @param offset
     * @param j
     */
    void LCP(
        NodeView &curr_node,
        int &offset,
        int &j);

    /**
     * @brief Check the end condition for the alignment.
     *
     * @param curr_data
     * @param curr_matrix  Matrix of the cell (M or MJumps)
     */
    void check_end_condition(
        Cell curr_data,
        Cell::Matrix curr_matrix);

    /**
     * @brief Exyend a given diagonal for a given vertex and perform the necessary
     * jumps.
     *
     * @param curr_node
     * @param prev_pos
     * @param prev_matrix
     */
    void extend_diagonal(
        NodeId curr_node_id,
        Cell::pos_t prev_pos,
        Cell::Matrix from_matrix);

    /**
     * @brief Initialize the partial backtrace, by finding the starting cell for backtrace.
     *
     */
    void init_partial_backtrace();

    /**
     * @brief Add matches to our backtracking vector.
     *
     * @param start_matches
     * @param end_matches
     */
    void add_matches(int start_matches, int end_matches);

    /**
     * @brief Add a mismatch to our backtracking vector.
     *
     */
    void add_mismatch();

    /**
     * @brief Add an insertion to our backtracking vector.
     *
     */
    void add_insertion();

    /**
     * @brief Add a deletion to our backtracking vector.
     *
     */
    void add_deletion();

    /**
     * @brief Perform a single step of the backtrace process.
     *
     * @param curr_cell
     * @param curr_matrix  Matrix of the current cell (updated with the one of the previous cell)
     */
    void one_backtrace_step(
        Cell &curr_cell,
        Cell::Matrix &curr_matrix);

    /**
     * @brief Backtrace the alignment from the end vertex to the start vertex.
     *
     */
    void backtrace();


    /**
     * @brief Internal function to print the graph in GFA format. This is used by
     * the public print_as_gfa function, which is the one that should be called by the user.
     *
     * @param gfa_output
     */
    void print_as_gfa_internal(std::ostream &gfa_output);

    /**
     * @brief Internal function to print the graph in dot (graphviz) format.
     *
     * @param out_stream Output stream to write the graph in dot format
     */
    void print_code_graphviz_internal(std::ostream &out_stream);

    // Handle reverse alignment
    NodeView get_node(NodeId id);
    bool has_out_nodes(NodeId id);

    int32_t _score = 0;

    Penalties _penalties;
    InternalPenalties _internal_penalties;

    std::unique_ptr<POAGraph> _poa_graph; // Partial order alignment graph for MSA

    bool _ends_free;
    bool _reversed_alignment = false;   // Read by the printing functions before the first alignment
    int  _start_column;
    int  _seq_ID = 0;
    NodeId _start_node;
    NodeId _end_node;
    int  _start_offset;
    Cell _start_pos;
    Cell::Matrix _start_matrix;
    int _start_score = 0;               // Internal score of the backtrace start cell

    std::unique_ptr<ScratchPad> _scratchpad;

    std::vector<ExtendState> _extend_stack;   // Reused by extend_diagonal

    std::unique_ptr<Scope> _scope;
    std::unique_ptr<BeyondScope> _beyond_scope;

    std::unique_ptr<VerticesData> _vertices_data;

    Heuristics _heuristics;

    Graph _graph;   // The graph to align to

    bool _is_msa;

    SequenceView _seq;

    Alignment _alignment;
};

}   // namespace theseus