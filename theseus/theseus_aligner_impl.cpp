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


#include <bit>
#include <cstring>
#include <string_view>
#include "theseus_aligner_impl.h"

namespace theseus {

TheseusAlignerImpl::TheseusAlignerImpl(const Penalties &penalties,
                                       const Heuristics &heuristics,
                                       Graph &&graph,
                                       int initial_weight,
                                       bool is_msa) :  _penalties(penalties),
                                                       _internal_penalties(penalties),
                                                       _heuristics(heuristics),
                                                       _graph(std::move(graph)),
                                                       _is_msa(is_msa),
                                                       _seq("", false) {
    // TODO: Gap-linear
    const auto n_scores = std::max({_internal_penalties.gapo() +_internal_penalties.dele(),
                                  _internal_penalties.gapo() +_internal_penalties.inse(),
                                  _internal_penalties.gapo2() +_internal_penalties.dele2(),
                                  _internal_penalties.gapo2() +_internal_penalties.inse2(),
                                  _internal_penalties.mism()}) + 1;

    // POA graph for MSA
    if (_is_msa) {
      _poa_graph = std::make_unique<POAGraph>();
      _poa_graph->create_initial_graph(_graph, initial_weight);
    }

    // Initialize aligner parameters and data structures
    _internal_penalties = InternalPenalties(penalties);
    _scope = std::make_unique<Scope>(n_scores);
    _beyond_scope = std::make_unique<BeyondScope>();
    constexpr int expected_nvertices = std::max(1024, 0); // TODO: Set the expected number of vertices
    _vertices_data = std::make_unique<VerticesData>(_internal_penalties, n_scores, expected_nvertices);
    _scratchpad = std::make_unique<ScratchPad>(-1024, 1024);
}

// Get the node/reversed node depending on the alignment configuration
NodeView TheseusAlignerImpl::get_node(NodeId id) {
  if (!_reversed_alignment) return _graph.node(id);
  else return _graph.node_rev(id);
}

bool TheseusAlignerImpl::has_out_nodes(NodeId id) {
  if (!_reversed_alignment) return !_graph.is_sink(id);
  else return !_graph.is_source(id);
}

void TheseusAlignerImpl::new_alignment(SequenceView seq,
                                       bool reverse_alignment,
                                       bool density_drop_active,
                                       bool lag_pruning_active) {
    _scope->new_alignment();
    _beyond_scope->new_alignment();
    _vertices_data->new_alignment();
    _seq = seq;
    _reversed_alignment = reverse_alignment;
    // Initialize scratchpad
    int max_diag = 0, v_n;
    NodeIdRange nodes = _graph.nodes();
    for (const auto &node : nodes) {
        v_n = _graph.node_size(node);
        max_diag = std::max(max_diag, v_n);
    }
    const int min_diag = -_seq.size();
    if (_scratchpad->max_diag() < max_diag ||
        _scratchpad->min_diag() > min_diag) {
        // TODO: Compute the max and min with a factor.
        _scratchpad = std::make_unique<ScratchPad>(min_diag, max_diag);
    }
    // Set data for first score
    _scope->new_score(_score);
    // Set initial alignment status
    _alignment.theseus_status = THESEUS_STATUS_OK;
    // Set heuristics
    _heuristics.new_alignment(std::max(_internal_penalties.inse(), _internal_penalties.dele()), _seq.size(), density_drop_active, lag_pruning_active);
    // TODO: Allow for different initial conditions. Now only global alignment.
    Cell init_condition;
    init_condition.offset    = 0;
    init_condition.vertex_id = _start_node;
    init_condition.diag      = _start_offset;
    init_condition.prev_pos  = -1;
    // Initial vertex data
    _beyond_scope->m_jumps_wf().push_back(init_condition);
    _vertices_data->activate_vertex(_start_node);
    _vertices_data->get_vertex_data(_start_node)._m_jumps_positions[0].push_back(0);
    // Alignment data
    _alignment.path.clear();
    _alignment.edit_op.clear();
    _alignment.score = 0;
}


// Process a given vertex with a given _score
void TheseusAlignerImpl::process_vertex(NodeId curr_node_id) {

  // Next
  int upper_bound = _graph.node_size(curr_node_id);
  const bool dual = _internal_penalties.is_dual();
  next_I<Cell::Matrix::I>(upper_bound, curr_node_id);
  _scratchpad->reset();
  if (dual) {
    next_I<Cell::Matrix::I2>(upper_bound, curr_node_id);
    _scratchpad->reset();
  }
  next_D<Cell::Matrix::D>(upper_bound, curr_node_id);
  _scratchpad->reset();
  if (dual) {
    next_D<Cell::Matrix::D2>(upper_bound, curr_node_id);
    _scratchpad->reset();
  }
  next_M(upper_bound, curr_node_id);
  _scratchpad->reset();
  // Extend. The vertex is fetched once and the diagonals that do not reach the
  // end of the vertex are extended in place (common case). Only those reaching
  // the end go through extend_diagonal, which performs the jumps.
  int v_pos = _vertices_data->get_id(curr_node_id);
  Scope::range cells_range = _scope->m_pos(_score)[v_pos];
  NodeView curr_node = get_node(curr_node_id);
  const int n = curr_node.sequence.size();
  for (Cell::pos_t idx = cells_range.start; idx < cells_range.end; ++idx) {
    // Local node extension
    Cell &cell = _beyond_scope->m_wf()[idx];  // Jumps never reallocate m_wf
    int j = cell.diag + cell.offset;
    LCP(curr_node, cell.offset, j);
    // If the diagonal reaches the end of the vertex, we need to jump to the
    // next vertices
    if (j == n) extend_diagonal(curr_node_id, idx, Cell::Matrix::M);
    else check_end_condition(cell, Cell::Matrix::M);
  }
}


void TheseusAlignerImpl::compute_new_wave() {
  // Update invalid segments. At score 0 the first jumps have just been created by
  // the initial extend (before this call), so the rem counters are already correct.
  if (_score > 0) _vertices_data->expand();
  _vertices_data->compact();
  // Process all active vertices
  int num_active_vertices = _vertices_data->num_active_vertices();
  NodeId curr_node_id;
  for (int l = 0; l < num_active_vertices; ++l) {
    curr_node_id = _vertices_data->get_vertex_id(l);
    process_vertex(curr_node_id);
  }
  _beyond_scope->update_positions();
}


Alignment TheseusAlignerImpl::align(
    std::string_view seq,
    // Seq-to-graph parameters
    int  start_node,
    int  start_offset,
    // MSA parameters
    int  weight,
    bool is_ends_free,
    // Common parameters
    bool reverse_alignment,
    bool density_drop_active,
    bool lag_pruning_active
  )
{
  // MSA_mode: Start position depends on the alignment direction)
  if (_is_msa) {
    if (!reverse_alignment) {
      _start_node = 0;
      _start_offset = 0;
      _end_node = 2;
    }
    else {
      _start_node = 2;
      _start_offset = 0;
      _end_node = 0;
    }
  }
  // Align a sequence to a reference graph starting at a given position
  else {
    _start_node   = start_node;
    _start_offset = start_offset;
  }
  // Set alignment parameters
  _ends_free = is_ends_free;
  // Initialize data for the new alignment
  new_alignment(Graph::SequenceView(seq, reverse_alignment), reverse_alignment, density_drop_active, lag_pruning_active);
  _score = 0;
  // _graph.print_code_graphviz();
  // Main alignment loop
  while (_alignment.theseus_status == THESEUS_STATUS_OK) {
    // Initial extend
    if (_score == 0) {
      extend_diagonal(_start_node, 0, Cell::Matrix::MJumps);
    }
    // Next wave + extend
    compute_new_wave();
    // Evaluate global heuristics
    _alignment.theseus_status = (_alignment.theseus_status == THESEUS_STATUS_ALG_COMPLETED) ?
                                 _alignment.theseus_status :
                                  _heuristics.check_global_heuristics(_score);
    // Update score
    _score = _score + 1;
    // Clear the corresponding waves and metadata from the scope
    _scope->new_score(_score);
    _vertices_data->new_score(_score);
  }
  _score -= 1;
  // Backtrace
  _seq_ID += 1;
  if (_alignment.theseus_status == THESEUS_STATUS_ALG_COMPLETED) {
    _start_score = _score;
    backtrace();
    if (_is_msa) {
      // Compute the end column of the alignment in the POA graph
      int end_column = _start_pos.offset + _start_pos.diag;
      _poa_graph->add_alignment_poa(_graph, _alignment, _seq, _seq_ID, weight, end_column);
    }
  }
  else {
    // Error messages
    if (_alignment.theseus_status == THESEUS_STATUS_END_UNREACHABLE) {
      // Choose starting position for backtrace (cell with maximum offset in M[s-_s_min]-M[s-_s_min - scope_size])
      // matrix on the last scope scores
      init_partial_backtrace();
      // Perform partial backtrace
      backtrace();
      // No drop in MSA mode
    }
  }
  return _alignment;
}

  // Sparsify M data
  void TheseusAlignerImpl::sparsify_M_data(Cell::CellVector & dense_wf,
                                           int offset_increase,
                                           int shift_factor,
                                           Scope::range cells_range,
                                           int m,
                                           int upper_bound)
  {
    Cell::pos_t len = cells_range.end - cells_range.start, new_col;
    Cell new_cell;
    // Sparsify the active diagonals
    for (int l = 0; l < len; ++l)
    {
      new_cell = dense_wf[cells_range.start + l];
      new_cell.diag += shift_factor;
      new_cell.offset += offset_increase;
      new_cell.from_matrix = Cell::Matrix::M;
      new_cell.prev_pos = cells_range.start + l;
      new_col = new_cell.offset + new_cell.diag; // d = j - i -> j = d + i
      // Check validity
      if (new_cell.offset <= m && new_col <= upper_bound)
      { // If in bounds
        // Branchless push_back
        auto &cell = _scratchpad->access_alloc(new_cell.diag);
        // If better offset
        const bool cmp = cell.offset < new_cell.offset;
        cell = (cmp) ? new_cell : cell;
      }
    }
  }

  // Sparsify jumps data
  void TheseusAlignerImpl::sparsify_jumps_data(Cell::CellVector & dense_wf,
                                               std::vector<Cell::pos_t> & jumps_positions,
                                               int offset_increase,
                                               int shift_factor,
                                               int m,
                                               int upper_bound,
                                               Cell::Matrix from_matrix)
  {
    int len = jumps_positions.size(), new_col, pos;
    Cell new_cell;
    // Sparsify the active diagonals
    for (int l = 0; l < len; ++l)
    {
      pos = jumps_positions[l];
      new_cell = dense_wf[pos];
      new_cell.prev_pos = pos;
      new_cell.from_matrix = from_matrix;
      new_cell.diag += shift_factor;
      new_cell.offset += offset_increase;
      new_col = new_cell.offset + new_cell.diag; // d = j - i -> j = d + i
      // Check validity
      if (new_cell.offset <= m && new_col <= upper_bound)
      { // If in bounds
        // Branchless push_back
        auto &cell = _scratchpad->access_alloc(new_cell.diag);
        // If better offset
        const bool cmp = cell.offset < new_cell.offset;
        cell = (cmp) ? new_cell : cell;
      }
    }
  }

  // Sparsify indel
  void TheseusAlignerImpl::sparsify_indel_data(Cell::CellVector & dense_wf,
                                               int offset_increase,
                                               int shift_factor,
                                               Scope::range cells_range,
                                               int m,
                                               int upper_bound)
  {
    Cell::pos_t len = cells_range.end - cells_range.start, new_col;
    Cell new_cell;
    // Sparsify the active diagonals
    for (int l = 0; l < len; ++l)
    {
      // Vertex_id and previous matrix are the same as before
      new_cell = dense_wf[cells_range.start + l];
      new_cell.diag += shift_factor;
      new_cell.offset += offset_increase;
      new_col = new_cell.offset + new_cell.diag; // d = j - i -> j = d + i
      // Check validity
      if (new_cell.offset <= m && new_col <= upper_bound)
      { // If in bounds
        // Branchless push_back
        auto &cell = _scratchpad->access_alloc(new_cell.diag);
        // If better offset
        const bool cmp = cell.offset < new_cell.offset;
        cell = (cmp) ? new_cell : cell;
      }
    }
  }

  // Data of the gap matrices, selected at compile time
  template <Cell::Matrix matrix>
  Cell::CellVector &TheseusAlignerImpl::gap_wf(int score) {
    if constexpr (matrix == Cell::Matrix::I)       return _scope->i_wf(score);
    else if constexpr (matrix == Cell::Matrix::I2) return _scope->i2_wf(score);
    else if constexpr (matrix == Cell::Matrix::D)  return _scope->d_wf(score);
    else                                           return _scope->d2_wf(score);
  }

  template <Cell::Matrix matrix>
  Scope::RangeVector &TheseusAlignerImpl::gap_pos(int score) {
    if constexpr (matrix == Cell::Matrix::I)       return _scope->i_pos(score);
    else if constexpr (matrix == Cell::Matrix::I2) return _scope->i2_pos(score);
    else if constexpr (matrix == Cell::Matrix::D)  return _scope->d_pos(score);
    else                                           return _scope->d2_pos(score);
  }

  template <Cell::Matrix matrix>
  Cell::CellVector &TheseusAlignerImpl::gap_jumps_wf() {
    if constexpr (matrix == Cell::Matrix::I) return _beyond_scope->i_jumps_wf();
    else                                     return _beyond_scope->i2_jumps_wf();
  }

  template <Cell::Matrix matrix>
  std::vector<std::vector<Cell::pos_t>> &TheseusAlignerImpl::gap_jumps_positions(NodeId node_id) {
    if constexpr (matrix == Cell::Matrix::I) return _vertices_data->get_vertex_data(node_id)._i_jumps_positions;
    else                                     return _vertices_data->get_vertex_data(node_id)._i2_jumps_positions;
  }

  // Compute next I (or I2) matrix
  template <Cell::Matrix matrix>
  void TheseusAlignerImpl::next_I(int upper_bound,
                                  NodeId curr_node_id)
  {
    // Penalties of the gap piece (I: first piece, I2: second piece)
    constexpr bool second = (matrix == Cell::Matrix::I2);
    const int gapo = second ? _internal_penalties.gapo2() : _internal_penalties.gapo();
    const int inse = second ? _internal_penalties.inse2() : _internal_penalties.inse();
    const int v_id = _vertices_data->get_id(curr_node_id);
    // Sparsify data (put it in the scratch pad)
    int pos_prev_M = _score - (gapo + inse),
        pos_prev_I = _score - inse;
    int pos_prev_M_scope = _vertices_data->get_pos(pos_prev_M),
        pos_prev_I_scope = _vertices_data->get_pos(pos_prev_I);
    // Come from an Insertion
    if (pos_prev_I >= 0) {
      if (gap_pos<matrix>(pos_prev_I).size() > v_id)
      {
        Scope::range cells_range = gap_pos<matrix>(pos_prev_I)[v_id];
        sparsify_indel_data(gap_wf<matrix>(pos_prev_I), 0, 1, cells_range, _seq.size(), upper_bound); // Sparsify I data
      };
      sparsify_jumps_data(
        gap_jumps_wf<matrix>(),
        gap_jumps_positions<matrix>(curr_node_id)[pos_prev_I_scope],
        0, 1, _seq.size(), upper_bound, second ? Cell::Matrix::I2Jumps : Cell::Matrix::IJumps);
    }
    // Come from M
    if (pos_prev_M >= 0) {
      if (_scope->m_pos(pos_prev_M).size() > v_id) {
        Scope::range cells_range = _scope->m_pos(pos_prev_M)[v_id];
        sparsify_M_data(_beyond_scope->m_wf(), 0, 1, cells_range, _seq.size(), upper_bound); // Sparsify M data
      }
      sparsify_jumps_data(_beyond_scope->m_jumps_wf(),
          _vertices_data->get_vertex_data(curr_node_id)._m_jumps_positions[pos_prev_M_scope],
          0, 1, _seq.size(), upper_bound, Cell::Matrix::MJumps);
    }
    // Densify data (store it in the big wavefront)
    Scope::range new_range;
    new_range.start = gap_wf<matrix>(_score).size();
    for (auto diag : _scratchpad->active_diags()) {
      if (_vertices_data->valid_diagonal<matrix>(curr_node_id, diag)
      && !_heuristics.check_local_heuristics((*_scratchpad)[diag].offset)) {
        gap_wf<matrix>(_score).push_back((*_scratchpad)[diag]);     // Store Cell
      }
    }
    new_range.end = gap_wf<matrix>(_score).size();
    gap_pos<matrix>(_score).push_back(new_range);
    // Check, store and invalidate new I jumps
    if (has_out_nodes(curr_node_id)) {
      NodeView curr_node = get_node(curr_node_id);
      check_and_store_jumps<matrix>(curr_node, gap_wf<matrix>(_score), new_range);
    }
}


// Compute next D (or D2) matrix
template <Cell::Matrix matrix>
void TheseusAlignerImpl::next_D(int upper_bound,
                                NodeId curr_node_id)
{
  // Penalties of the gap piece (D: first piece, D2: second piece)
  constexpr bool second = (matrix == Cell::Matrix::D2);
  const int gapo = second ? _internal_penalties.gapo2() : _internal_penalties.gapo();
  const int dele = second ? _internal_penalties.dele2() : _internal_penalties.dele();
  const int v_id = _vertices_data->get_id(curr_node_id);
  // Sparsify data (put it in the scratch pad)
  int pos_prev_M = _score - (gapo + dele),
      pos_prev_D = _score - dele,
      pos_prev_M_scope = _vertices_data->get_pos(pos_prev_M);
  // Come from a Deletion
  if (pos_prev_D >= 0 && gap_pos<matrix>(pos_prev_D).size() > v_id)
  {
    Scope::range cells_range = gap_pos<matrix>(pos_prev_D)[v_id];
    sparsify_indel_data(gap_wf<matrix>(pos_prev_D), 1, -1, cells_range, _seq.size(), upper_bound);
  }
  // Come from M
  if (pos_prev_M >= 0) {
    if (_scope->m_pos(pos_prev_M).size() > v_id)
    {
      Scope::range cells_range = _scope->m_pos(pos_prev_M)[v_id];
      sparsify_M_data(_beyond_scope->m_wf(), 1, -1, cells_range, _seq.size(), upper_bound); // Sparsify M data
    }
    sparsify_jumps_data(_beyond_scope->m_jumps_wf(),
        _vertices_data->get_vertex_data(curr_node_id)._m_jumps_positions[pos_prev_M_scope],
        1, -1, _seq.size(), upper_bound, Cell::Matrix::MJumps);
  }
  // Densify data (store it in the big wavefront)
  Scope::range new_range;
  new_range.start = gap_wf<matrix>(_score).size();
  for (auto diag : _scratchpad->active_diags()) {
    if (_vertices_data->valid_diagonal<matrix>(curr_node_id, diag)
    && !_heuristics.check_local_heuristics((*_scratchpad)[diag].offset)) {
      gap_wf<matrix>(_score).push_back((*_scratchpad)[diag]); // Store Cell
    }
  }
  new_range.end = gap_wf<matrix>(_score).size();
  gap_pos<matrix>(_score).push_back(new_range);
}


// Compute next M matrix
void TheseusAlignerImpl::next_M(int upper_bound,
                                NodeId curr_node_id) {
  // Sparsify data (put it in the scratch pad)
  int pos_prev_M = _score - _internal_penalties.mism(),
      pos_prev_D = _score,
      pos_prev_I = _score,
      pos_prev_M_scope = _vertices_data->get_pos(pos_prev_M);
  // Come from a Deletion
  if (_scope->d_pos(pos_prev_D).size() > _vertices_data->get_id(curr_node_id))  {
    Scope::range cells_range = _scope->d_pos(pos_prev_D)[_vertices_data->get_id(curr_node_id)];
    sparsify_indel_data(_scope->d_wf(pos_prev_D), 0, 0, cells_range, _seq.size(), upper_bound);
  }
  // Come from an Insertion
  if (_scope->i_pos(pos_prev_I).size() > _vertices_data->get_id(curr_node_id))  {
    Scope::range cells_range = _scope->i_pos(pos_prev_I)[_vertices_data->get_id(curr_node_id)];
    sparsify_indel_data(_scope->i_wf(pos_prev_I), 0, 0, cells_range, _seq.size(), upper_bound);
  }
  // Come from the second piece of the gaps (dual affine)
  if (_internal_penalties.is_dual()) {
    if (_scope->d2_pos(pos_prev_D).size() > _vertices_data->get_id(curr_node_id))  {
      Scope::range cells_range = _scope->d2_pos(pos_prev_D)[_vertices_data->get_id(curr_node_id)];
      sparsify_indel_data(_scope->d2_wf(pos_prev_D), 0, 0, cells_range, _seq.size(), upper_bound);
    }
    if (_scope->i2_pos(pos_prev_I).size() > _vertices_data->get_id(curr_node_id))  {
      Scope::range cells_range = _scope->i2_pos(pos_prev_I)[_vertices_data->get_id(curr_node_id)];
      sparsify_indel_data(_scope->i2_wf(pos_prev_I), 0, 0, cells_range, _seq.size(), upper_bound);
    }
  }
  // Come from M
  if (pos_prev_M >= 0) {
    if (_scope->m_pos(pos_prev_M).size() > _vertices_data->get_id(curr_node_id))  {
      Scope::range cells_range = _scope->m_pos(pos_prev_M)[_vertices_data->get_id(curr_node_id)];
      sparsify_M_data(_beyond_scope->m_wf(), 1, 0, cells_range,  _seq.size(), upper_bound);
    }
    sparsify_jumps_data(_beyond_scope->m_jumps_wf(),
        _vertices_data->get_vertex_data(curr_node_id)._m_jumps_positions[pos_prev_M_scope],
        1, 0, _seq.size(), upper_bound, Cell::Matrix::MJumps);
  }
  // Densify data (store it in the big wavefront)
  Scope::range new_range;
  new_range.start = _beyond_scope->m_wf().size();
  for (auto diag : _scratchpad->active_diags()) {
    if (_vertices_data->valid_diagonal<Cell::Matrix::M>(curr_node_id, diag)
    && !_heuristics.check_local_heuristics((*_scratchpad)[diag].offset)) {
      _beyond_scope->m_wf().push_back((*_scratchpad)[diag]);     // Store Cell
    }
  }
  new_range.end = _beyond_scope->m_wf().size();
  _scope->m_pos(_score).push_back(new_range);
}


// Store the jump in neighbours. "jump_cell" is the cell that reached the end of
// its vertex (a copy, as the underlying vectors might be reallocated)
void TheseusAlignerImpl::store_M_jump(NodeView curr_node,
                                      Cell jump_cell) {
  // Invalidate the jumping diagonal
  _vertices_data->invalidate_m_jump(_vertices_data->get_id(jump_cell.vertex_id), jump_cell.diag);
  int pos_score = _vertices_data->get_pos(_score);
  // The jump cell keeps from_matrix and prev_pos (the origin of the gap) for the backtrace
  jump_cell.diag = -jump_cell.offset;
  // For each neighbour, store the jump and metadata
  for (auto out_node_id : curr_node.out_nodes) {
    jump_cell.vertex_id = out_node_id;
    _vertices_data->activate_vertex(jump_cell.vertex_id);
    // Store jump and metadata
    bool valid_diag = _vertices_data->valid_diagonal<Cell::Matrix::M>(jump_cell.vertex_id, jump_cell.diag);
    // Extend only if it has not yet been visited
    if (valid_diag) {
      int pos_new_cell = _beyond_scope->m_jumps_wf().size();
      _beyond_scope->m_jumps_wf().push_back(jump_cell);
      _vertices_data->get_vertex_data(jump_cell.vertex_id)._m_jumps_positions[pos_score].push_back(pos_new_cell);
      // The new M cell (start of the neighbour) must be extended and checked
      extend_diagonal(out_node_id, pos_new_cell, Cell::Matrix::MJumps);
    }
  }
}


// Store the jump in a destination vertex. If the destination is empty, its
// neighbours are added to the pending jumps
template <Cell::Matrix matrix>
void TheseusAlignerImpl::store_I_jump_in_vertex(
    NodeId dest_id,
    Cell::pos_t prev_pos,
    Cell::Matrix from_matrix,
    Cell &jump_cell,
    int pos_score,
    std::vector<PendingIJump> &pending)
{
  constexpr Cell::Matrix jumps_matrix = (matrix == Cell::Matrix::I2) ? Cell::Matrix::I2Jumps : Cell::Matrix::IJumps;
  jump_cell.vertex_id   = dest_id;
  jump_cell.prev_pos    = prev_pos;
  jump_cell.from_matrix = from_matrix;
  _vertices_data->activate_vertex(dest_id);
  bool valid_diag = _vertices_data->valid_diagonal<matrix>(dest_id, jump_cell.diag);
  // Extend only if it has not yet been visited
  if (valid_diag) {
    int pos_new_cell = gap_jumps_wf<matrix>().size();
    gap_jumps_wf<matrix>().push_back(jump_cell);
    gap_jumps_positions<matrix>(dest_id)[pos_score].push_back(pos_new_cell);
    // If the destination vertex is empty, jump again
    NodeView dest_node = get_node(dest_id);
    if (dest_node.sequence.empty() && has_out_nodes(dest_id)) {
      if constexpr (matrix == Cell::Matrix::I2) {
        _vertices_data->invalidate_i2_jump(_vertices_data->get_id(dest_id), jump_cell.diag);
      }
      else {
        _vertices_data->invalidate_i_jump(_vertices_data->get_id(dest_id), jump_cell.diag);
      }
      for (auto next_node_id : dest_node.out_nodes) {
        pending.push_back({next_node_id, pos_new_cell, jumps_matrix});
      }
    }
  }
}


// Store the jump in neighbours (insertion that continues in the neighbours)
template <Cell::Matrix matrix>
void TheseusAlignerImpl::store_I_jump(
    NodeView curr_node,
    Cell jump_cell)
{
  // Invalidate the jumping diagonal
  if constexpr (matrix == Cell::Matrix::I2) {
    _vertices_data->invalidate_i2_jump(_vertices_data->get_id(jump_cell.vertex_id), jump_cell.diag);
  }
  else {
    _vertices_data->invalidate_i_jump(_vertices_data->get_id(jump_cell.vertex_id), jump_cell.diag);
  }
  int pos_score = _vertices_data->get_pos(_score);
  jump_cell.diag = -jump_cell.offset;
  // If a destination vertex is empty the insertion keeps going through it, so its
  // neighbours become pending jumps (explicit stack instead of recursion). The
  // cell stored in the empty vertex is the origin of the next one. The vector only
  // allocates memory if there are empty vertices.
  std::vector<PendingIJump> pending;
  // For each neighbour, store the jump and metadata
  const Cell::pos_t origin_pos = jump_cell.prev_pos;
  const Cell::Matrix origin_matrix = jump_cell.from_matrix;
  for (auto out_node_id : curr_node.out_nodes) {
    store_I_jump_in_vertex<matrix>(out_node_id, origin_pos, origin_matrix, jump_cell, pos_score, pending);
  }
  // Pending jumps through empty vertices
  while (!pending.empty()) {
    PendingIJump curr_jump = pending.back();
    pending.pop_back();
    store_I_jump_in_vertex<matrix>(curr_jump.node_id, curr_jump.prev_pos, curr_jump.from_matrix, jump_cell, pos_score, pending);
  }
}


// Check and store I (or I2) jumps (that is, those diagonals that have reached the last column of a vertex)
template <Cell::Matrix matrix>
void TheseusAlignerImpl::check_and_store_jumps(
    NodeView curr_node,
    Cell::CellVector &curr_wavefront,
    Scope::range cell_range)
{
  Cell::pos_t len = cell_range.end - cell_range.start, n = curr_node.sequence.size();
  // Check all diagonals in the current wavefront
  for (int l = 0; l < len; ++l) {
    const Cell cell = curr_wavefront[cell_range.start + l];
    // Reached the end of the vertex
    if (cell.diag + cell.offset == n && cell.offset <= (int)_seq.size()) {
      store_M_jump(curr_node, cell);
      store_I_jump<matrix>(curr_node, cell);
    }
  }
}


// Compute the Longest Common Prefix between two given sequences
void TheseusAlignerImpl::LCP(
    NodeView &curr_node,
    int &offset,
    int &j)
{
  // Find LCP
  int len_seq_1 = _seq.size();
  int len_seq_2 = curr_node.sequence.size();
  // Compare 8 characters at a time: the first different byte of the XOR of the
  // two blocks is the first mismatch. In a reversed view, the next 8 characters
  // are the 8 bytes before the current position, in reverse order.
  if constexpr (std::endian::native == std::endian::little) {
    const char *seq_1 = _seq.data(), *seq_2 = curr_node.sequence.data();
    const bool reversed = _seq.is_reversed();
    if (reversed == curr_node.sequence.is_reversed()) {
      while (offset + 8 <= len_seq_1 && j + 8 <= len_seq_2) {
        uint64_t block_1, block_2;
        std::memcpy(&block_1, reversed ? seq_1 + len_seq_1 - 8 - offset : seq_1 + offset, 8);
        std::memcpy(&block_2, reversed ? seq_2 + len_seq_2 - 8 - j : seq_2 + j, 8);
        const uint64_t diff = block_1 ^ block_2;
        if (diff != 0) {
          const int equal = (reversed ? std::countl_zero(diff) : std::countr_zero(diff)) / 8;
          offset += equal;
          j += equal;
          return;
        }
        offset += 8;
        j += 8;
      }
    }
  }
  // Remaining characters, one at a time
  while (offset < len_seq_1 && j < len_seq_2 && _seq[offset] == curr_node.sequence[j]) {
    offset = offset + 1;   // Update the f.r. of this diagonal
    j = j + 1;
  }
}


// TODO: Implement different end conditions as Global, Semi-Global...
void TheseusAlignerImpl::check_end_condition(
    Cell curr_data,
    Cell::Matrix curr_matrix)
{
  if (curr_data.offset == (int)_seq.size()) {
    _alignment.theseus_status = THESEUS_STATUS_ALG_COMPLETED;
    _start_pos    = curr_data;
    _start_matrix = curr_matrix;
  }
}


// Extend a particular diagonal
void TheseusAlignerImpl::extend_diagonal(
    NodeId       init_node_id,
    Cell::pos_t  init_pos,
    Cell::Matrix init_matrix)
{
  // Use a explicit stack to avoid recursion and stack overflow. The stack is a
  // member so its memory is reused (allocating it for every diagonal was the
  // most expensive part of the extend). extend_diagonal is never reentered.
  std::vector<ExtendState> &extend_stack = _extend_stack;
  // Push the initial state onto the stack
  extend_stack.push_back({init_node_id, init_pos, init_matrix});
  // Process the stack until it's empty
  while (!extend_stack.empty()) {
    auto [curr_node_id, curr_pos, curr_from_matrix] = extend_stack.back();
    extend_stack.pop_back();

    // Extend the current diagonal
    NodeView curr_node_view = get_node(curr_node_id);
    Cell* curr_cell_ptr = (curr_from_matrix == Cell::Matrix::M) ? &_beyond_scope->m_wf()[curr_pos] :
                          &_beyond_scope->m_jumps_wf()[curr_pos];
    int j = curr_cell_ptr->diag + curr_cell_ptr->offset;
    LCP(curr_node_view, curr_cell_ptr->offset, j);
    // End condition
    check_end_condition(*curr_cell_ptr, curr_from_matrix);

    // Jump to neighbours if the end of the current node is reached
    if (j == (int)curr_node_view.sequence.size() && curr_cell_ptr->offset <= (int)_seq.size() && has_out_nodes(curr_node_id)) {
      // Invalidate the jumping diagonal
      _vertices_data->invalidate_m_jump(_vertices_data->get_id(curr_cell_ptr->vertex_id), curr_cell_ptr->diag);
      // Compute data of the cell in the next vertex
      int pos_score = _vertices_data->get_pos(_score);
      Cell new_cell = *curr_cell_ptr;
      new_cell.from_matrix = curr_from_matrix;
      new_cell.prev_pos    = curr_pos;
      new_cell.diag        = -curr_cell_ptr->offset;
      // For the corresponding neighbour, store the jump and metadata
      for (auto out_node_id : curr_node_view.out_nodes) {
        new_cell.vertex_id = out_node_id;
        _vertices_data->activate_vertex(new_cell.vertex_id);
        // Store jump and metadata
        bool valid_diag = _vertices_data->valid_diagonal<Cell::Matrix::M>(new_cell.vertex_id, new_cell.diag);
        // Extend only if it has not yet been visited
        if (valid_diag) {
          int pos_new_cell = _beyond_scope->m_jumps_wf().size();
          _beyond_scope->m_jumps_wf().push_back(new_cell);
          _vertices_data->get_vertex_data(new_cell.vertex_id)._m_jumps_positions[pos_score].push_back(pos_new_cell);
          // Push the next state onto the stack for the next neighbour
          extend_stack.push_back({out_node_id, pos_new_cell, Cell::Matrix::MJumps});
        }
      }
    }
  }
}


// Initialize the partial backtrace, by finding the starting cell for backtrace.
// We find this cell by checking all the active cells in scores (s-_s_min, ..., s-_s_min-_scope_size)
void TheseusAlignerImpl::init_partial_backtrace() {
  Cell best_cell;
  best_cell.offset = -1;
  // Iterate through the valid score range
  int start_score = _score - _heuristics.s_min() - _scope->size();
  start_score = std::max(0, start_score);
  // Iterate in the M structure for score s TODO: I_jumps?
  for (int s = start_score; s <= _score - _heuristics.s_min(); ++s) {
    // Check M wavefront
    int start_pos_M = (s == 0) ? 0 : _beyond_scope->m_wf_pos(s-1);
    int end_pos_M   = _beyond_scope->m_wf_pos(s);   // Exclusive
    for (int pos = start_pos_M; pos < end_pos_M; ++pos) {
      Cell curr_cell = _beyond_scope->m_wf()[pos];
      if (curr_cell.offset > best_cell.offset) {
        best_cell = curr_cell;
        _start_score = s;
        _start_matrix = Cell::Matrix::M;
      }
    }
    // Check M_jumps wavefront
    int start_pos_M_jumps = (s == 0) ? 0 : _beyond_scope->m_jumps_wf_pos(s-1);
    int end_pos_M_jumps   = _beyond_scope->m_jumps_wf_pos(s);   // Exclusive
    for (int pos = start_pos_M_jumps; pos < end_pos_M_jumps; ++pos) {
      Cell curr_cell = _beyond_scope->m_jumps_wf()[pos];
      if (curr_cell.offset > best_cell.offset) {
        best_cell = curr_cell;
        _start_score = s;
        _start_matrix = Cell::Matrix::MJumps;
      }
    }
  }
  _start_pos = best_cell;
}


// Add matches to our backtracking vector
void TheseusAlignerImpl::add_matches(
    int start_matches,
    int end_matches)
{
  int size = end_matches - start_matches;
  for (int k = 0; k < size; ++k) {
    _alignment.edit_op.push_back('M');
  }
}


// Add a mismatch to our backtracking vector
void TheseusAlignerImpl::add_mismatch()
{
  _alignment.edit_op.push_back('X');
}


// Add an insertion to our backtracking vector
void TheseusAlignerImpl::add_insertion()
{
  _alignment.edit_op.push_back('I');
}


// Add a deletion to our backtracking vector
void TheseusAlignerImpl::add_deletion()
{
  _alignment.edit_op.push_back('D');
}


void TheseusAlignerImpl::one_backtrace_step(
    Cell &curr_cell,
    Cell::Matrix &curr_matrix)
{
  Cell prev_cell;
  if (curr_cell.from_matrix == Cell::Matrix::M) prev_cell = _beyond_scope->m_wf()[curr_cell.prev_pos];
  else if (curr_cell.from_matrix == Cell::Matrix::MJumps) prev_cell = _beyond_scope->m_jumps_wf()[curr_cell.prev_pos];
  else if (curr_cell.from_matrix == Cell::Matrix::I2Jumps) prev_cell = _beyond_scope->i2_jumps_wf()[curr_cell.prev_pos];
  else prev_cell = _beyond_scope->i_jumps_wf()[curr_cell.prev_pos];
  int num_indels;
  // We jump if the current matrix is MJumps, IJumps or I2Jumps. Comparing the
  // vertices is not enough, as a vertex can jump to itself (self-loop)
  bool is_jump = (curr_matrix == Cell::Matrix::MJumps || curr_matrix == Cell::Matrix::IJumps || curr_matrix == Cell::Matrix::I2Jumps);
  // We are inside the same vertex
  if (!is_jump) { // Still in the same vertex
    // Mismatch
    if (curr_cell.diag == prev_cell.diag) {
      if (curr_cell.offset > prev_cell.offset) {    // Consider 0 length vertices
        add_matches(prev_cell.offset + 1, curr_cell.offset);
        add_mismatch();
      }
    }
    else {
      // Deletion
      if (curr_cell.diag < prev_cell.diag) {
        num_indels = prev_cell.diag - curr_cell.diag;
        add_matches(prev_cell.offset + num_indels, curr_cell.offset);
        for (int l = 0; l < num_indels; ++l) add_deletion();
      }
      // Insertion
      else {
        num_indels = curr_cell.diag - prev_cell.diag;
        add_matches(prev_cell.offset, curr_cell.offset);
        for (int l = 0; l < num_indels; ++l) add_insertion();
      }
    }
  }
  // Jump to another vertex
  else {
    add_matches(prev_cell.offset, curr_cell.offset);                          // Add the necessary matches
    _alignment.path.push_back(prev_cell.vertex_id);                           // Add the new vertex to the path
    int col_in_prev_v = prev_cell.diag + prev_cell.offset;
    int num_insertions = _graph.node_size(prev_cell.vertex_id) - col_in_prev_v;
    for (int l = 0; l < num_insertions; ++l) add_insertion();                 // Add the necessary insertions
  }
  // Update current matrix and cell
  curr_matrix = curr_cell.from_matrix;
  curr_cell = prev_cell;
}


// Main function of the backtracking process
void TheseusAlignerImpl::backtrace()
{
  Cell curr_pos = _start_pos;
  _alignment.start_offset = _start_offset;
  _alignment.score = _internal_penalties.to_original_score(_start_score, _start_pos.offset);
  _alignment.end_offset = curr_pos.diag + curr_pos.offset; // Vertex offset = j
  _alignment.path.push_back(curr_pos.vertex_id);
  Cell::Matrix curr_matrix = _start_matrix;
  // Main backtrace loop
  while (curr_pos.prev_pos != -1)
  {
    one_backtrace_step(curr_pos, curr_matrix);
  }
  // Add the matches until the beginning of the sequence
  add_matches(0, curr_pos.offset);
  // Reverse path and edit operations
  if (!_reversed_alignment) {
    std::reverse(_alignment.edit_op.begin(), _alignment.edit_op.end());
    std::reverse(_alignment.path.begin(), _alignment.path.end());
  }
}


// Output functions
/**
 * @brief Visualize the graph in Graphviz format.
 *
 * @param G
 */
void TheseusAlignerImpl::print_code_graphviz_internal(std::ostream &out_stream)
{
    // TODO: What if we wanted to print reversed?
    _reversed_alignment = false;
    out_stream << "digraph G {" << std::endl;
    // Print nodes
    for (NodeId id : _graph.nodes())
    {
      NodeView node = get_node(id);
      out_stream << id << " [label=\"";
      for (size_t j = 0; j < node.sequence.size(); ++j)
      {
          out_stream << node.sequence[j];
      }
      out_stream << "\"]" << std::endl;
    }
    // Print edges
    for (NodeId id : _graph.nodes())
    {
      auto node = get_node(id);
      for (NodeId out_id : node.out_nodes)
      {
        out_stream << id << "->" << out_id << std::endl;
      }
    }
    out_stream << "}" << std::endl;
}
void TheseusAlignerImpl::print_code_graphviz(std::ostream &out_stream) {
  print_code_graphviz_internal(out_stream);
}

// Print as GFA
void TheseusAlignerImpl::print_as_gfa_internal(
  std::ostream &gfa_output)
{
    // Print all nodes as Segments
    for (NodeId id : _graph.nodes())
    {
      NodeView node = get_node(id);
      gfa_output << "S\t" << id << "\t" << node.sequence << "\n";
    }
    // Print all edges as Links
    for (NodeId id : _graph.nodes())
    {
      NodeView node = get_node(id);
      // Go through all incoming vertices (with this you cover all possible edges,
      // since the graph is directed)
      for (const auto in_node : node.in_nodes)
      {
        gfa_output << "L\t" << in_node << "\t+\t"
          << id << "\t+\t" << "0M\n";
      }
    }
}
// TODO: Using node names?
void TheseusAlignerImpl::print_as_gfa(std::ostream &gfa_output)
{
  print_as_gfa_internal(gfa_output);
}

// Print as msa (can only call from TheseusMSA)
void TheseusAlignerImpl::print_as_msa(std::ostream &out_stream) {
  _poa_graph->poa_to_fasta(_seq_ID, out_stream);
}


// Find and return the consensus sequence (can only call from TheseusMSA)
std::string TheseusAlignerImpl::heaviest_bundle_consensus() {
  return _poa_graph->poa_to_consensus();
}

// Find and return the consensus sequence using weighted majority voting
void TheseusAlignerImpl::majority_voting_consensus(std::vector<int> &consensus_weights,
                                                   std::string &consensus_sequence,
                                                   std::string &consensus_sequence_gapped) {
  _poa_graph->poa_to_consensus_weighted_majority_voting(_seq_ID, consensus_weights, consensus_sequence, consensus_sequence_gapped);
}

// Print as GAF
void TheseusAlignerImpl::print_as_gaf(
    theseus::Alignment &alignment,
    std::ostream &out_stream,
    std::string seq_name,
    std::unordered_map<NodeId, std::string> &node_names) {

  // Field 1: Query name
  out_stream << seq_name;

  // Field 2: Query length
  out_stream << "\t" << _seq.size();

  // Field 3: Query start
  out_stream << "\t" << 0;

  // Field 4: Query end
  out_stream << "\t" << _seq.size();

  // Field 5: Strand
  out_stream << "\t" << "+"; // TODO: Support reverse strand

  // Field 6: Alignment path
  out_stream << "\t";
  for (size_t l = 0; l < alignment.path.size(); ++l) {
    out_stream << ">" << node_names.at(alignment.path[l]); // TODO: Support orientation
  }

  // Field 7: Target length
  int target_length = 0;
  for (size_t l = 0; l < alignment.path.size(); ++l) {
    int node_length = node_names.at(alignment.path[l]).size();
    target_length += node_length;
  }
  out_stream << "\t" << target_length;

  // Field 8: Target start
  out_stream << "\t" << _start_offset;

  // Field 9: Target end
  out_stream << "\t" << alignment.end_offset;

  // Field 10: Number of matching bases
  int num_matches = 0;
  for (size_t l = 0; l < alignment.edit_op.size(); ++l) {
    if (alignment.edit_op[l] == 'M') {
      num_matches += 1;
    }
  }
  out_stream << "\t" << num_matches;

  // Field 11: Alignment block length
  out_stream << "\t" << alignment.edit_op.size();

  // Field 12: Mapping quality
  out_stream << "\t" << 255; // TODO: Compute mapping quality

  // Optional fields
  out_stream << "\t" << "cg:Z:"; // CIGAR string
  std::string cigar = "";
  int count = 1;
  for (size_t l = 1; l < alignment.edit_op.size(); ++l) {
    if (alignment.edit_op[l] == alignment.edit_op[l - 1]) {
      count += 1;
    }
    else {
      cigar += std::to_string(count) + alignment.edit_op[l - 1];
      count = 1;
    }
  }
  if (alignment.edit_op.size() > 0) {
    cigar += std::to_string(count) + alignment.edit_op.back();
  }
  out_stream << cigar << "\n";
}

} // namespace theseus