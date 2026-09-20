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

#include "../doctest.h"

#include "../../include/theseus/penalties.h"

TEST_CASE("Check penalty members") {
    SUBCASE("Linear") {
        theseus::Penalties penalties(0, 1, 1);
        CHECK(penalties.type() == theseus::Penalties::Linear);
        CHECK(penalties.match() == 0);
        CHECK(penalties.mism() == 1);
        CHECK(penalties.gape() == 1);
        CHECK(penalties.gapo() == 0);
    }

    SUBCASE("Affine") {
        theseus::Penalties penalties(0, 1, 2, 1);
        CHECK(penalties.type() == theseus::Penalties::Affine);
        CHECK(penalties.match() == 0);
        CHECK(penalties.mism() == 1);
        CHECK(penalties.gapo() == 2);
        CHECK(penalties.gape() == 1);
        CHECK(penalties.gapo2() == 0);
    }

    SUBCASE("DualAffine") {
        theseus::Penalties penalties(0, 1, 2, 2, 5, 1);
        CHECK(penalties.type() == theseus::Penalties::DualAffine);
        CHECK(penalties.match() == 0);
        CHECK(penalties.mism() == 1);
        CHECK(penalties.gapo() == 2);
        CHECK(penalties.gape() == 2);
        CHECK(penalties.gapo2() == 5);
        CHECK(penalties.gape2() == 1);
    }
}