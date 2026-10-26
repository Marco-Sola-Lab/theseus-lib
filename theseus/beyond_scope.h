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

#include <vector>

#include "cell.h"
#include "vector.h"

/**
 * Class containing the necessary wavefronts to perform backtrace. These wavefronts
 * have to be stored in memory until the end of the alignment.
 *
 */

namespace theseus {

class BeyondScope {
public:
    /**
     * @brief Construct a new Beyond Scope object
     *
     */
    BeyondScope() {
        constexpr int expected_ncells = 4096;

        _m_wf.realloc(expected_ncells);
        _m_jumps_wf.realloc(expected_ncells);
        _i_jumps_wf.realloc(expected_ncells);
        _i2_jumps_wf.realloc(expected_ncells);

        _m_wf.set_realloc_policy(dense_wf_realloc_policy);
        _m_jumps_wf.set_realloc_policy(dense_wf_realloc_policy);
        _i_jumps_wf.set_realloc_policy(dense_wf_realloc_policy);
        _i2_jumps_wf.set_realloc_policy(dense_wf_realloc_policy);
    }

    /**
     * @brief Reinitialize the beyond the scope object each time that a new
     * alignment is called.
     *
     */
    void new_alignment() {
        _m_wf.clear();
        _m_jumps_wf.clear();
        _i_jumps_wf.clear();
        _i2_jumps_wf.clear();

        _m_wf_pos.clear();
        _m_jumps_wf_pos.clear();
        _i_jumps_wf_pos.clear();
        _i2_jumps_wf_pos.clear();
    }

    /**
     * @brief Update the positions vectors. This means, compute the size of the
     * Cell vectors (M, M_jumps, I_jumps and I2_jumps) and push this value back to
     * the associated vectors.
     *
     */
    void update_positions() {
        _m_wf_pos.push_back(_m_wf.size());
        _m_jumps_wf_pos.push_back(_m_jumps_wf.size());
        _i_jumps_wf_pos.push_back(_i_jumps_wf.size());
        _i2_jumps_wf_pos.push_back(_i2_jumps_wf.size());
    }


    /**
     * @brief Access the i_jumps wavefront
     *
     * @return Cell::Wavefront&
     */
    Cell::CellVector &i_jumps_wf() {
        return _i_jumps_wf;
    }

    /**
     * @brief Access the i2_jumps wavefront
     *
     * @return Cell::Wavefront&
     */
    Cell::CellVector &i2_jumps_wf() {
        return _i2_jumps_wf;
    }

    /**
     * @brief Access the m_jumps wavefront
     *
     * @return Cell::Wavefront&
     */
    Cell::CellVector &m_jumps_wf() {
        return _m_jumps_wf;
    }

    /**
     * @brief Access the m wavefront
     *
     * @return Cell::Wavefront&
     */
    Cell::CellVector &m_wf() {
        return _m_wf;
    }

    /**
     * @brief Access the m wavefront positions vector
     *
     * @return Cell::Wavefront&
     */
    int m_wf_pos(int s) {
        return _m_wf_pos[s];
    }

    /**
     * @brief Access the m jumps wavefront positions vector
     *
     */
    int m_jumps_wf_pos(int s) {
        return _m_jumps_wf_pos[s];
    }

    /**
     * @brief Access the i jumps wavefront positions vector
     *
     */
    int i_jumps_wf_pos(int s) {
        return _i_jumps_wf_pos[s];
    }

    /**
     * @brief Access the i2 jumps wavefront positions vector
     *
     */
    int i2_jumps_wf_pos(int s) {
        return _i2_jumps_wf_pos[s];
    }


private:
    static constexpr std::ptrdiff_t dense_wf_realloc_policy([[maybe_unused]] std::ptrdiff_t capacity,
                                                            std::ptrdiff_t required_size) {
        return required_size * 2;
    };

    Cell::CellVector _m_wf;        // M structure backtrace wavefront
    Cell::CellVector _m_jumps_wf;  // M Jumps structure backtrace wavefront
    Cell::CellVector _i_jumps_wf;  // I Jumps structure backtrace wavefront
    Cell::CellVector _i2_jumps_wf; // I2 Jumps structure backtrace wavefront

    // Vectors to know the span of cells associated to each score
    std::vector<int> _m_wf_pos;
    std::vector<int> _m_jumps_wf_pos;
    std::vector<int> _i_jumps_wf_pos;
    std::vector<int> _i2_jumps_wf_pos;
};

} // namespace theseus