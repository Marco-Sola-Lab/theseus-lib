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

#include "cell.h"
#include <vector>

/**
 * The scope class manages the temporary wavefront data used during alignment.
 * That is, it stores the wavefronts and position vectors for each score in a
 * circular queue.
 *
 * The M and jumps wavefronts are normally kept for the whole alignment (for the
 * backtrace) in the BeyondScope. In score-only mode they are not needed after
 * leaving the scope, so they are stored here too (see m_wf and the jumps
 * wavefronts below).
 *
 */

namespace theseus {

class Scope {
public:
    using pos_t = int64_t;

    /**
     * @brief Struct representing a range of values
     *
     */
    struct range {
        pos_t start;
        pos_t end;
    };

    // TODO: Prefer this?
    using RangeVector = Vector<range, true>;
    // using JumpsPos = ManualCapacityVector<int32_t>;

    /**
     * @brief Construct a new Scope object
     *
     * @param nscores Number of scores in the scope
     */
    Scope(int nscores) {
        _squeue.realloc(nscores);

        for (int i = 0; i < nscores; i++) {
            constexpr int init_capacity = 1024;
            ScoreData sd(init_capacity);
            _squeue.push_back(std::move(sd));
        }
    }

    /**
     * @brief Restore data for a new alignment.
     *
     */
    void new_alignment() {
        for (int i = 0; i < _squeue.size(); ++i) {
            _squeue[i].resize(0);
        }
    }

    /**
     * @brief Reinitialize data for a new score. As the scope works as a circular
     * queue, the previously stored data in a position of the scope is not relevant
     * anymore and its contents should be reset to be reused again.
     *
     */
    void new_score(int score) {
        _squeue[score%_squeue.size()].resize(0);
    }

    /**
     * @brief Get the size of the scope.
     *
     * @return int
     */
    int size() {
        return _squeue.size();
    }


    /**
     * @brief Get the data from the wavefront I of score "score".
     *
     * @param score
     * @return Cell::CellVector&
     */
    Cell::CellVector &i_wf(int score) {
        return _squeue[score%_squeue.size()]._i_wf;
    }

    /**
     * @brief Get the data from the wavefront D of score "score".
     *
     * @param score
     * @return Cell::CellVector&
     */
    Cell::CellVector &d_wf(int score) {
        return _squeue[score%_squeue.size()]._d_wf;
    }

    /**
     * @brief Get the data from the wavefront I2 of score "score".
     *
     * @param score
     * @return Cell::CellVector&
     */
    Cell::CellVector &i2_wf(int score) {
        return _squeue[score%_squeue.size()]._i2_wf;
    }

    /**
     * @brief Get the data from the wavefront D2 of score "score".
     *
     * @param score
     * @return Cell::CellVector&
     */
    Cell::CellVector &d2_wf(int score) {
        return _squeue[score%_squeue.size()]._d2_wf;
    }

    /**
     * @brief Get the cells of the M wavefronts of score "score" (score-only
     * mode only).
     *
     * @param score
     * @return Cell::CellVector&
     */
    Cell::CellVector &m_wf(int score) {
        return _squeue[score%_squeue.size()]._m_wf;
    }

    /**
     * @brief Get the cells of the M jumps of score "score" (score-only mode
     * only).
     *
     * @param score
     * @return Cell::CellVector&
     */
    Cell::CellVector &m_jumps_wf(int score) {
        return _squeue[score%_squeue.size()]._m_jumps_wf;
    }

    /**
     * @brief Get the cells of the I jumps of score "score" (score-only mode
     * only).
     *
     * @param score
     * @return Cell::CellVector&
     */
    Cell::CellVector &i_jumps_wf(int score) {
        return _squeue[score%_squeue.size()]._i_jumps_wf;
    }

    /**
     * @brief Get the cells of the I2 jumps of score "score" (score-only mode
     * only).
     *
     * @param score
     * @return Cell::CellVector&
     */
    Cell::CellVector &i2_jumps_wf(int score) {
        return _squeue[score%_squeue.size()]._i2_jumps_wf;
    }

    /**
     * @brief Get the data from the vector of M positions at score "score".
     *
     * @param score
     * @return RangeVector&
     */
    RangeVector &m_pos(int score) {
        return _squeue[score%_squeue.size()]._m_pos;
    }

    /**
     * @brief Get the data from the vector of I positions at score "score".
     *
     * @param score
     * @return RangeVector&
     */
    RangeVector &i_pos(int score) {
        return _squeue[score%_squeue.size()]._i_pos;
    }

    /**
     * @brief Get the data from the vector of I2 positions at score "score".
     *
     * @param score
     * @return RangeVector&
     */
    RangeVector &i2_pos(int score) {
        return _squeue[score%_squeue.size()]._i2_pos;
    }

    /**
     * @brief Get the data from the vector of D positions at score "score".
     *
     * @param score
     * @return RangeVector&
     */
    RangeVector &d_pos(int score) {
        return _squeue[score%_squeue.size()]._d_pos;
    }

    /**
     * @brief Get the data from the vector of D2 positions at score "score".
     *
     * @param score
     * @return RangeVector&
     */
    RangeVector &d2_pos(int score) {
        return _squeue[score%_squeue.size()]._d2_pos;
    }

private:
    struct ScoreData {
        static constexpr std::ptrdiff_t realloc_policy([[maybe_unused]] std::ptrdiff_t capacity,
                                                       std::ptrdiff_t required_size) {
            return required_size * 2;
        };

        Cell::CellVector _i_wf;
        Cell::CellVector _d_wf;

        Cell::CellVector _i2_wf;
        Cell::CellVector _d2_wf;

        // M and jumps cells (score-only mode only, otherwise they are kept in
        // the BeyondScope). They start empty and grow on demand.
        Cell::CellVector _m_wf;
        Cell::CellVector _m_jumps_wf;
        Cell::CellVector _i_jumps_wf;
        Cell::CellVector _i2_jumps_wf;

        RangeVector _m_pos;

        RangeVector _i_pos;
        RangeVector _i2_pos;

        RangeVector _d_pos;
        RangeVector _d2_pos;

        ScoreData(int capacity) {
            _i_wf.realloc(capacity);
            _d_wf.realloc(capacity);
            _i2_wf.realloc(capacity);
            _d2_wf.realloc(capacity);

            _m_pos.realloc(capacity);
            _i_pos.realloc(capacity);
            _i2_pos.realloc(capacity);
            _d_pos.realloc(capacity);
            _d2_pos.realloc(capacity);

            _i_wf.set_realloc_policy(realloc_policy);
            _d_wf.set_realloc_policy(realloc_policy);
            _i2_wf.set_realloc_policy(realloc_policy);
            _d2_wf.set_realloc_policy(realloc_policy);

            _m_wf.set_realloc_policy(realloc_policy);
            _m_jumps_wf.set_realloc_policy(realloc_policy);
            _i_jumps_wf.set_realloc_policy(realloc_policy);
            _i2_jumps_wf.set_realloc_policy(realloc_policy);

            _m_pos.set_realloc_policy(realloc_policy);
            _i_pos.set_realloc_policy(realloc_policy);
            _i2_pos.set_realloc_policy(realloc_policy);
            _d_pos.set_realloc_policy(realloc_policy);
            _d2_pos.set_realloc_policy(realloc_policy);
        }

        void resize(int new_size) {
            _i_wf.resize(new_size);
            _d_wf.resize(new_size);
            _i2_wf.resize(new_size);
            _d2_wf.resize(new_size);

            _m_wf.resize(new_size);
            _m_jumps_wf.resize(new_size);
            _i_jumps_wf.resize(new_size);
            _i2_jumps_wf.resize(new_size);

            _m_pos.resize(new_size);
            _i_pos.resize(new_size);
            _i2_pos.resize(new_size);
            _d_pos.resize(new_size);
            _d2_pos.resize(new_size);
        }
    };

    Vector<ScoreData> _squeue;
};

} // namespace theseus