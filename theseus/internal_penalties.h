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

#include <stdexcept>

#include "theseus/penalties.h"

/**
 * A class derived from Penalties that contains the internal penalties used
 * during alignment.
 *
 */

namespace theseus {

class InternalPenalties {
public:
    using penalty_t = Penalties::penalty_t;

    InternalPenalties(Penalties penalties) : _original_match(penalties.match()) {
        // Transform penalties to have match = 0
        if (penalties.match() != 0) {
            // What happens with negative or mixed penalties? TODO:
            // We subtract "match" once per consumed QUERY character (M, X and D).
            // The query length is fixed, so every path in the graph is shifted by
            // the same constant. Subtracting it per graph character (I) would
            // depend on the length of the chosen path in the graph, potentilly
            // yielding different, non-equivalent, optimal paths.
            _match    = 0;
            _mismatch = penalties.mism() - penalties.match();
            _gapo     = penalties.gapo();
            _inse     = penalties.gape();
            _dele     = penalties.gape() - penalties.match();
        }
        else {
            _match    = penalties.match();
            _mismatch = penalties.mism();
            _gapo     = penalties.gapo();
            _inse     = penalties.gape();
            _dele     = penalties.gape();
        }

        // Second piece of the gap function
        _dual = penalties.type() == Penalties::Type::DualAffine;
        if (_dual) {
            _gapo2 = penalties.gapo2();
            _inse2 = penalties.gape2();
            _dele2 = penalties.gape2() - penalties.match();
        }
        else {
            _gapo2 = _gapo;
            _inse2 = _inse;
            _dele2 = _dele;
        }

        // We ask the user for minimization penalties
        // TODO: Accept maximization penalties and internally transform them to
        // minimization penalties
        if (penalties.match() >= penalties.mism()) {
            throw std::invalid_argument("The match penalty must be less than the mismatch penalty");
        }
        else if (penalties.match() >= penalties.gapo() + penalties.gape()) {
            throw std::invalid_argument("The gap open + gap extension penalty must be greater than the match penalty.");
        }
        else if (penalties.match() >= penalties.gape()) {
            throw std::invalid_argument("The gap extension penalty must be greater than the match penalty.");
        }
        else if (penalties.gapo() < penalties.gape()) {
            throw std::invalid_argument("The gap open penalty must be greater than or equal to the gap extension penalty.");
        }
        else if (penalties.gapo() + penalties.gape() <= 0) {
            throw std::invalid_argument("The sum of the gap open and gap extension penalties must be greater than 0.");
        }
        else if (penalties.gape() <= 0) {
            throw std::invalid_argument("The gap extension penalty must be greater than 0.");
        }
        // Dual affine: the second piece must be the "long gaps" piece (higher
        // opening and lower extension).
        if (_dual) {
            if (penalties.gape2() <= penalties.match()) {
                throw std::invalid_argument("The second gap extension penalty must be greater than the match penalty.");
            }
            else if (penalties.gape2() <= 0) {
                throw std::invalid_argument("The second gap extension penalty must be greater than 0.");
            }
            else if (penalties.gape2() > penalties.gape()) {
                throw std::invalid_argument("The second gap extension penalty must be less than or equal to the first one.");
            }
            else if (penalties.gapo2() < penalties.gapo()) {
                throw std::invalid_argument("The second gap open penalty must be greater than or equal to the first one.");
            }
        }
    }

    /**
     * Transform an internal score back to the score under the user penalties.
     * Internal penalties subtract "match" once per consumed query character in
     * every alignment path (also with dual affine gaps), so:
     *   original = internal + match * consumed_chars
     * Valid for full alignments (consumed_chars = sequence length) and partial ones.
     *
     * @param internal_score Score under the internal penalties.
     * @param consumed_chars Number of query characters consumed.
     * @return The score under the user penalties.
     */
    penalty_t to_original_score(penalty_t internal_score, int consumed_chars) const {
        return internal_score + _original_match * consumed_chars;
    }

    /**
     * Whether the gap function has two pieces (dual affine).
     *
     * @return True if the penalties are dual affine.
     */
    bool is_dual() const { return _dual; }

    /**
     * Get the match score.
     *
     * @return The match score.
     */
    penalty_t match() const { return _match; }

    /**
     * Get the mismatch score.
     *
     * @return The mismatch score.
     */
    penalty_t mism() const { return _mismatch; }

    /**
     * Get the gap open penalty if the gap type is affine or dual affine.
     *
     * @return The gap open penalty if the gap type is affine or dual affine.
     */
    penalty_t gapo() const { return _gapo; }

    /**
     * Get the insertion extension penalty.
     *
     * @return The insertion extension penalty.
     */
    penalty_t inse() const { return _inse; }

    /**
     * Get the deletion extension penalty.
     *
     * @return The deletion extension penalty.
     */
    penalty_t dele() const { return _dele; }

    /**
     * Get the second gap open penalty if the gap type is dual affine.
     *
     * @return The gap open penalty if the gap type is dual affine.
     */
    penalty_t gapo2() const { return _gapo2; }

    /**
     * Get the second insertion extension penalty if the gap type is dual affine.
     *
     * @return The insertion extension penalty if the gap type is dual affine.
     */
    penalty_t inse2() const { return _inse2; }

    /**
     * Get the second deletion extension penalty if the gap type is dual affine.
     *
     * @return The deletion extension penalty if the gap type is dual affine.
     */
    penalty_t dele2() const { return _dele2; }

private:
    penalty_t _original_match; // Match score of the user penalties

    // Diagonal penalties
    penalty_t _match;
    penalty_t _mismatch;

    // Affine gap penalties
    penalty_t _gapo;
    penalty_t _inse;
    penalty_t _dele;

    // Dual affine gap penalties
    penalty_t _gapo2;
    penalty_t _inse2;
    penalty_t _dele2;
    bool _dual;
};

} // namespace theseus