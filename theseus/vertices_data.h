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

#include <algorithm>
#include <cstdint>
#include <vector>
// #include <type_traits>

#include "cell.h"
#include "theseus/graph.h"
#include "internal_penalties.h"


namespace theseus {

using NodeId = Graph::NodeId;

/**
 * @brief Class vertices data. It stores the data related to the active vertices,
 * their indexes, the invalid diagonals and some jumping data.
 *
 */
class VerticesData {
public:
    using pos_t = int64_t;

    struct Segment {
        int32_t start_d;   // Start diagonal of the segment (inclusive)
        int32_t end_d;     // End diagonal of the segment (inclusive)
    };

    struct InvalidData {
        Segment seg;        // Segment of invalid diagonals
        int32_t rem_up;     // Remaining scores to grow the segment one diagonal
                            // up
        int32_t rem_down;   // Remaining scores to grow the segment one diagonal
                            // down
    };


    /**
     * @brief Vertex dat structure. It contains:
     * - Invalid data for the M, I (I2) and D (D2) matrices.
     * - Positions of the jumps in the MJ and IJ (I2) matrices.
     */
    struct VertexData {
        NodeId vertex_id;

        std::vector<InvalidData> _m_invalid;

        std::vector<InvalidData> _i_invalid;
        std::vector<InvalidData> _i2_invalid;

        std::vector<InvalidData> _d_invalid;
        std::vector<InvalidData> _d2_invalid;

        // Scope with the positions of M jumps in the scope previous waves
        std::vector<std::vector<pos_t>> _m_jumps_positions;

        // Scope with the positions of I jumps in the scope previous waves
        std::vector<std::vector<pos_t>> _i_jumps_positions;

        // Scope with the positions of I2s jumps in the scope previous waves
        std::vector<std::vector<pos_t>> _i2_jumps_positions;
    };

    int _nscores;

    /**
     * @brief Construct a new Vertices Data object
     *
     * @param penalties Set of used penalties.
     * @param nscores   Number of scores of the scope.
     * @param nexpected_vertices    Number of expected vertices.
     */
    VerticesData(const InternalPenalties &penalties, int nscores, int nexpected_vertices) :
        _nscores(nscores), _penalties(penalties) {
        _active_vertices.reserve(nexpected_vertices);
        _vertex_to_idx.reserve(nexpected_vertices);
    }

    /**
     * @brief Clear the previously stored values in the scope to write the
     * newly relevant data.
     *
     * @param score  Current score
     */
    void new_score(int score) {
        int len = _active_vertices.size();
        int pos_curr_score = get_pos(score);

        for (int l = 0; l < len; ++l) {
            auto &vdata = _active_vertices[l];

            // Clear the jumps (they work as a scope)
            vdata._m_jumps_positions[pos_curr_score].clear();
            vdata._i_jumps_positions[pos_curr_score].clear();
            vdata._i2_jumps_positions[pos_curr_score].clear();
        }
    }

    /**
     * @brief Reinitialize the vertices data object each time that a new
     * alignment is called.
     *
     */
    void new_alignment() {
        _active_vertices.clear();
        _vertex_to_idx.clear();
    }

    /**
     * @brief Get the id of a given vertex.
     *
     * @param vtx   Vertex id
     * @return Cell::vertex_t  Index of the vertex in the active vertices
     */
    int get_id(int vtx) {
        return _vertex_to_idx[vtx];
    }

    /**
     * @brief Get the vertex id of a given vertex.
     *
     * @param idx   Index of the vertex in the active vertices
     * @return Cell::vertex_t  Vertex id
     */
    Cell::vertex_t get_vertex_id(int idx) {
        return _active_vertices[idx].vertex_id;
    }

    /**
     * @brief Get the vertex data of a given vertex.
     *
     * @param vtx  Vertex id
     * @return VertexData&  Data of the vertex
     */
    VertexData &get_vertex_data(int vtx) {
        return _active_vertices[_vertex_to_idx[vtx]];
    }

    /**
     * @brief Get the position in the scope of a given score.
     *
     * @param score
     * @return int  Position in the scope (it ranges [0, _nscores-1])
     */
    int get_pos(int score) {
        return score%_nscores;
    }

    /**
     * @brief Get the number of scores.
     *
     * @return int Number of scores
     */
    int get_n_scores() {
        return _nscores;
    }

    // FUNCTIONS
    /**
     * @brief Compact a set of ordered invalid objects to avoid redundant information.
     *
     */
    void compact_invalid_vector(std::vector<InvalidData> &invalid_v,
                                int default_rem_up,
                                int default_rem_down) {

        if (invalid_v.size() == 0) {
            return;
        }

        // Order the set of segments
        std::sort(invalid_v.begin(),
                  invalid_v.end(),
                  [](const InvalidData &s1, const InvalidData &s2) {
                      return s1.seg.start_d < s2.seg.start_d;
                  });

        // Iterate through the loop
        int k = 0;
        for (long unsigned int l = 1; l < invalid_v.size(); ++l) {
            // You should compact them
            if (invalid_v[k].seg.end_d + 1 >= invalid_v[l].seg.start_d) {
                // Update remaining up scores. This must be done BEFORE updating the
                // end of the segment, as the comparison uses the old end of k.
                if (invalid_v[l].seg.end_d > invalid_v[k].seg.end_d) {
                    invalid_v[k].rem_up = std::min(invalid_v[l].rem_up,
                                                   invalid_v[k].rem_up +
                                                       (invalid_v[l].seg.end_d -
                                                        invalid_v[k].seg.end_d) *
                                                           default_rem_up);
                }
                else {
                    invalid_v[k].rem_up = std::min(invalid_v[k].rem_up,
                                                   invalid_v[l].rem_up +
                                                       (invalid_v[k].seg.end_d -
                                                        invalid_v[l].seg.end_d) *
                                                           default_rem_up);
                }

                // Update remaining down scores
                invalid_v[k].rem_down = std::min(invalid_v[k].rem_down,
                                                 invalid_v[l].rem_down +
                                                     (invalid_v[l].seg.start_d -
                                                      invalid_v[k].seg.start_d) *
                                                         default_rem_down);

                // Update segment
                invalid_v[k].seg.end_d = std::max(invalid_v[l].seg.end_d,
                                                  invalid_v[k].seg.end_d);
            }
            else {
                k += 1;
                invalid_v[k] = invalid_v[l];
            }
        }
        invalid_v.resize(k + 1);
    }

    /**
     * @brief Expand a set of invalid objects. This means reducing the counters
     * of remaining scores to grow and growing the segments if those counters get
     * to 0.
     *
     */
    void expand_invalid_vector(std::vector<InvalidData> &invalid_v,
                               int default_rem_up,
                               int default_rem_down) {

        for (long unsigned int l = 0; l < invalid_v.size(); ++l) {
            invalid_v[l].rem_down -= 1;
            invalid_v[l].rem_up -= 1;

            if (invalid_v[l].rem_up == 0) {
                invalid_v[l].rem_up = default_rem_up;
                invalid_v[l].seg.end_d += 1;
            }
            if (invalid_v[l].rem_down == 0) {
                invalid_v[l].rem_down = default_rem_down;
                invalid_v[l].seg.start_d -= 1;
            }
        }
    }

    /**
     * @brief Expand all invalid objects.
     *
     */
    void expand() {
        for (long unsigned int l = 0; l < _active_vertices.size(); ++l) {
            auto &vdata = _active_vertices[l];
            expand_invalid_vector(vdata._m_invalid, _penalties.inse(), _penalties.dele());
            expand_invalid_vector(vdata._i_invalid, _penalties.inse(), _penalties.dele());
            expand_invalid_vector(vdata._d_invalid, _penalties.inse(), _penalties.dele());

            // The growth rates of I2 and D2 are those of their own piece in the
            // direction of the gap (up for I2, down for D2). The other direction
            // is reached through a gap of the first piece.
            expand_invalid_vector(vdata._i2_invalid, _penalties.inse2(), _penalties.dele());
            expand_invalid_vector(vdata._d2_invalid, _penalties.inse(), _penalties.dele2());
        }
    }

    /**
     * @brief Compact all invalid objects.
     *
     */
    void compact() {
        for (long unsigned int l = 0; l < _active_vertices.size(); ++l) {
            compact_invalid_vector(_active_vertices[l]._m_invalid,
                                   _penalties.inse(),
                                   _penalties.dele());
            compact_invalid_vector(_active_vertices[l]._i_invalid,
                                   _penalties.inse(),
                                   _penalties.dele());
            compact_invalid_vector(_active_vertices[l]._d_invalid,
                                   _penalties.inse(),
                                   _penalties.dele());

            compact_invalid_vector(_active_vertices[l]._i2_invalid,
                                   _penalties.inse2(),
                                   _penalties.dele());
            compact_invalid_vector(_active_vertices[l]._d2_invalid,
                                   _penalties.inse(),
                                   _penalties.dele2());
        }
    }

    /**
     * @brief Invalidate a diagonal "diag" in vertex located at index "idx" and
     * for matrix I. This happens when a jump is performed in the I matrix.
     *
     * @param idx
     * @param diag
     */
    void invalidate_i_jump(int idx, int diag) {
        VertexData &vdata = _active_vertices[idx];
        const int open_i1 = _penalties.gapo() + _penalties.inse();
        const int open_d1 = _penalties.gapo() + _penalties.dele();
        const int open_i2 = _penalties.gapo2() + _penalties.inse2();
        const int open_d2 = _penalties.gapo2() + _penalties.dele2();

        // New invalid in M and I
        add_invalid(vdata._m_invalid, diag, diag, _penalties.inse(), open_d1);
        add_invalid(vdata._i_invalid, diag, diag, _penalties.inse(),
                    open_d1 + _penalties.dele() + open_i1);
        // New invalid in D (initially empty)
        add_invalid(vdata._d_invalid, diag, diag - 1, _penalties.inse() + open_d1, open_d1);

        if (_penalties.is_dual()) {
            // The I jump does not dominate the I2/D2 states (they extend cheaper),
            // so these segments start empty. Reaching I2 (D2) from the jump means
            // closing the gap and opening a new gap of the second piece.
            add_invalid(vdata._i2_invalid, diag + 1, diag, open_i2, open_d1 + open_i2);
            add_invalid(vdata._d2_invalid, diag, diag - 1, _penalties.inse() + open_d2, open_d2);
        }
    }

    /**
     * @brief Invalidate a diagonal "diag" in vertex located at index "idx" and
     * for matrix I2. This happens when a jump is performed in the I2 matrix.
     * Dual affine only. The I2 state dominates M, I and I2 (its extension is the
     * cheapest one), but not the D states.
     *
     * @param idx
     * @param diag
     */
    void invalidate_i2_jump(int idx, int diag) {
        VertexData &vdata = _active_vertices[idx];
        const int open_i1 = _penalties.gapo() + _penalties.inse();
        const int open_d1 = _penalties.gapo() + _penalties.dele();
        const int open_i2 = _penalties.gapo2() + _penalties.inse2();
        const int open_d2 = _penalties.gapo2() + _penalties.dele2();

        // New invalid in M, I and I2
        add_invalid(vdata._m_invalid, diag, diag, _penalties.inse2(), open_d1);
        add_invalid(vdata._i_invalid, diag, diag, open_i1,
                    open_d1 + _penalties.dele() + open_i1);
        add_invalid(vdata._i2_invalid, diag, diag, _penalties.inse2(),
                    open_d1 + _penalties.dele() + open_i2);
        // New invalid in D and D2 (initially empty)
        add_invalid(vdata._d_invalid, diag, diag - 1, _penalties.inse2() + open_d1, open_d1);
        add_invalid(vdata._d2_invalid, diag, diag - 1, _penalties.inse2() + open_d2, open_d2);
    }

    /**
     * @brief Invalidate a diagonal "diag" in vertex located at index "idx" and
     * for matrix M. This happens when a jump is performed in the M matrix.
     *
     * @param idx
     * @param diag
     */
    void invalidate_m_jump(int idx, int diag) {
        VertexData &vdata = _active_vertices[idx];
        const int open_i1 = _penalties.gapo() + _penalties.inse();
        const int open_d1 = _penalties.gapo() + _penalties.dele();
        const int open_i2 = _penalties.gapo2() + _penalties.inse2();
        const int open_d2 = _penalties.gapo2() + _penalties.dele2();

        // New invalid in M
        add_invalid(vdata._m_invalid, diag, diag, open_i1, open_d1);
        // New invalid in I and D (initially empty)
        add_invalid(vdata._i_invalid, diag + 1, diag, open_i1, open_d1 + open_i1);
        add_invalid(vdata._d_invalid, diag, diag - 1, open_i1 + open_d1, open_d1);

        if (_penalties.is_dual()) {
            // I2 and D2 start empty. They are reached by a gap of the second piece
            // (preceded by a gap of the first piece if it goes in the other direction)
            add_invalid(vdata._i2_invalid, diag + 1, diag, open_i2, open_d1 + open_i2);
            add_invalid(vdata._d2_invalid, diag, diag - 1, open_i1 + open_d2, open_d2);
        }
    }

    /**
     * @brief Check if a given diagonal "diag" is valid in the vertex located at
     * position "idx". This function allows to validate the M, I and D matrices.
     *
     * @tparam matrix
     * @param idx
     * @param diag
     * @return bool
     */
    template <Cell::Matrix matrix>
    bool valid_diagonal(int vtx, int diag) {
        std::vector<InvalidData> &invalid =
        [this, vtx]() -> std::vector<InvalidData>& {
            VertexData &vdata = _active_vertices[get_id(vtx)];
            if constexpr (matrix == Cell::Matrix::M) {
                return vdata._m_invalid;
            }
            else if constexpr (matrix == Cell::Matrix::I) {
                return vdata._i_invalid;
            }
            else if constexpr (matrix == Cell::Matrix::D) {
                return vdata._d_invalid;
            }
            else if constexpr (matrix == Cell::Matrix::I2) {
                return vdata._i2_invalid;
            }
            else if constexpr (matrix == Cell::Matrix::D2) {
                return vdata._d2_invalid;
            }
            else {
                static_assert([]{ return false; }(), "Unsupported matrix type");
            }
        }();

        for (size_t l = 0; l < invalid.size(); ++l) {
            if (invalid[l].seg.start_d <= diag &&
                diag <= invalid[l].seg.end_d) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Return the number of active vertices.
     *
     * @return int
     */
    size_t num_active_vertices() {
        return _active_vertices.size();
    }

    /**
     * @brief Activate a new vertex "vtx" in the active vertices list. Nothing is
     * done if the vertex is already active.
     *
     * @param vtx
     */
    void activate_vertex(NodeId vtx) {
        if (_vertex_to_idx.size() <= vtx) {
            _vertex_to_idx.resize(2*vtx + 1, -1);
        }
        if (_vertex_to_idx[vtx] == -1) {
            // Add the vertex to the active vertices
            _active_vertices.push_back(VertexData());
            _active_vertices[_active_vertices.size() - 1].vertex_id = vtx;
            _active_vertices[_active_vertices.size() - 1]._i_jumps_positions.resize(_nscores);
            _active_vertices[_active_vertices.size() - 1]._i2_jumps_positions.resize(_nscores);
            _active_vertices[_active_vertices.size() - 1]._m_jumps_positions.resize(_nscores);

            // Determine the vertex id
            _vertex_to_idx[vtx] = _active_vertices.size() - 1;
        }
    }

private:
    /**
     * @brief Add a new invalid segment [start_d, end_d] to an invalid vector.
     *
     */
    static void add_invalid(std::vector<InvalidData> &invalid_v,
                            int start_d,
                            int end_d,
                            int rem_up,
                            int rem_down) {
        InvalidData new_invalid;
        new_invalid.seg.start_d = start_d;
        new_invalid.seg.end_d = end_d;
        new_invalid.rem_up = rem_up;
        new_invalid.rem_down = rem_down;
        invalid_v.push_back(new_invalid);
    }

    const InternalPenalties &_penalties;

    std::vector<VertexData> _active_vertices;

    std::vector<int> _vertex_to_idx;
};

}   // namespace theseus