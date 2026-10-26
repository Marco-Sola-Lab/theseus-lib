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


#include "theseus/theseus_aligner.h"

#include "theseus_aligner_impl.h"

namespace theseus {

TheseusAligner::TheseusAligner(const Penalties &penalties,
                               const Heuristics &heuristics,
                               Graph &&graph)
{
    aligner_impl_ = std::make_unique<TheseusAlignerImpl>(penalties, heuristics, std::move(graph), 1, false);
}

TheseusAligner::TheseusAligner(const Penalties &penalties,
                               const Heuristics &heuristics,
                               Graph &graph)
{
    Graph graph_copy = graph; // Make a copy of the graph to ensure that the aligner has its own instance
    aligner_impl_ = std::make_unique<TheseusAlignerImpl>(penalties, heuristics, std::move(graph_copy), 1, false);
}

TheseusAligner::~TheseusAligner() {}

void TheseusAligner::print_alignment_as_gaf(
                theseus::Alignment &alignment,
                std::ostream &out_stream,
                std::string seq_name,
                std::unordered_map<NodeId, std::string> &node_names) {

    aligner_impl_->print_as_gaf(alignment, out_stream, seq_name, node_names);
}

/**
 * @brief Main alignment function for the Theseus aligner.
 *
 * @param seq
 * @param start_node
 * @param start_offset
 * @param density_drop_active
 * @param lag_pruning_active
 * @return Alignment
 */
Alignment TheseusAligner::align(
    std::string_view seq,
    NodeId &start_node,
    int start_offset,
    bool density_drop_active,
    bool lag_pruning_active) {

    // Arguments: weight, is_ends_free, reverse_alignment, density_drop_active, lag_pruning_active
    return aligner_impl_->align(seq, start_node, start_offset, 1, false, false, density_drop_active, lag_pruning_active);
}

void TheseusAligner::set_alignment_mode(AlignmentMode mode) {
    aligner_impl_->set_alignment_mode(mode);
}

} // namespace theseus