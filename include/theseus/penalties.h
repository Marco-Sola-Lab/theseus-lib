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

#include <span>
#include <string>
#include <vector>

/**
 * Class containing the penalties for the alignment algorithm. The objective
 * function is to minimize the alignment score.
 *
 */

namespace theseus {

class Penalties {
public:
    using penalty_t = int;

    enum Type {
        Linear,
        Affine,
        DualAffine
    };

    /**
     * Create a Linear-Gaps penalty object.
     *
     * @param match The match score.
     * @param mismatch The mismatch score.
     * @param gape The gap extension penalty.
     */
    Penalties(penalty_t match,
              penalty_t mismatch,
              penalty_t gape);

    /**
     * Create an Affine-Gap penalty object.
     *
     * @param match The match score.
     * @param mismatch The mismatch score.
     * @param gapo The gap open penalty.
     * @param gape The gap extension penalty.
     */
    Penalties(penalty_t match,
              penalty_t mismatch,
              penalty_t gapo,
              penalty_t gape);

    /**
     * Create a Dual Affine-Gap penalty object.
     *
     * @param match The match score.
     * @param mismatch The mismatch score.
     * @param gapo The first gap open penalty.
     * @param gape The first gap extension penalty.
     * @param gapo2 The second gap open penalty.
     * @param gape2 The second gap extension penalty.
     */
    Penalties(penalty_t match,
              penalty_t mismatch,
              penalty_t gapo,
              penalty_t gape,
              penalty_t gapo2,
              penalty_t gape2);

    /**
     * Get the gap type.
     *
     * @return The gap type.
     */
    Type type() const { return type_; }

    /**
     * Get the match score.
     *
     * @return The match score.
     */
    penalty_t match() const { return match_; }

    /**
     * Get the mismatch score.
     *
     * @return The mismatch score.
     */
    penalty_t mism() const { return mismatch_; }

    /**
     * Get the gap open penalty if the gap type is affine or dual affine.
     * Otherwise, return 0.
     *
     * @return The gap open penalty if the gap type is affine or dual affine.
     * Otherwise, return 0.
     */
    penalty_t gapo() const { return gapo_; }

    /**
     * Get the gap extension penalty.
     *
     * @return The gap extension penalty.
     */
    penalty_t gape() const { return gape_; }

    /**
     * Get the second gap open penalty if the gap type is dual affine.
     * Otherwise, return 0.
     *
     * @return The gap open penalty if the gap type is dual affine. Otherwise,
     * return 0.
     */
    penalty_t gapo2() const { return gapo2_; }

    /**
     * Get the second gap extension penalty if the gap type is dual affine.
     * Otherwise, return 0.
     *
     * @return The gap extension penalty if the gap type is dual affine.
     * Otherwise, return 0.
     */
    penalty_t gape2() const { return gape2_; }


protected:
    Type type_;

    penalty_t match_; // Match score.
    penalty_t mismatch_; // Mismatch score.

    penalty_t gapo_; // Gap open.
    penalty_t gape_; // Gap extension.

    penalty_t gapo2_; // Gap open for dual gap-affine.
    penalty_t gape2_; // Gap extension for dual gap-affine.
};

} // namespace theseus
