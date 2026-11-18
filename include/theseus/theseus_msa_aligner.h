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
#include <istream>
#include <unordered_map>

#include "theseus/penalties.h"
#include "theseus/alignment.h"
#include "theseus/heuristics.h"
#include "theseus/graph.h"

/**
 * Multiple Sequence Aligner (MSA) based on POA graphs. Internally uses the
 * TheseusAlignerImpl class.
 *
 */

namespace theseus
{

    class TheseusAlignerImpl; // Forward declaration of the implementation class.

    class TheseusMSA
    {
    public:
        /**
         * Initial constructor
         *
         * @param penalties    User defined alignment penalties
         * @param heuristics   User defined heuristics
         * @param seq          Sequence to initialize the graph
         * @param initial_weight    Initial weight for the first added sequence
         */
        TheseusMSA(
            const Penalties &penalties,
            const Heuristics &heuristics,
            std::string_view seq,
            int initial_weight);

        /**
         * Class destructor
         *
         */
        ~TheseusMSA();

        /**
         * Add a new sequence to the POA graph, representing the MSA so far.
         * The sequence is aligned until the whole query is consumed (see
         * EndCondition::FullQuery) and it is always linked to the sink.
         *
         * @param seq Sequence to add to the MSA
         * @param weight Weight of the sequence
         * @param lag_pruning_active Whether to use the lag pruning heuristic
         */
        Alignment align(std::string_view seq,
                        int  weight = 1,
                        bool lag_pruning_active = false);

        /**
         * Add a new sequence to the POA graph with the options of this sequence:
         * lag pruning, reverse alignment, end condition and density drop (see
         * MSAOptions). If the alignment is dropped by the density drop heuristic,
         * the partial alignment is added to the graph and the status of the
         * returned alignment is THESEUS_STATUS_END_UNREACHABLE.
         *
         * @param seq Sequence to add to the MSA
         * @param weight Weight of the sequence
         * @param options Options of the alignment of this sequence
         */
        Alignment align(std::string_view seq,
                        int  weight,
                        const MSAOptions &options);

        /**
         * @brief Print the current POA graph as a GFA file.
         *
         */
        void print_as_gfa(std::ostream &out_stream);


        /**
         * @brief Print the current POA graph in MSA format.
         *
         */
        void print_as_msa(std::ostream &out_stream);


        /**
         * @brief Return consensus sequence.
         *
         */
        std::string heaviest_bundle_consensus();

        /**
         * @brief Compute the weighted majority voting consensus sequence.
         *
         */
        void majority_voting_consensus(std::vector<int> &consensus_weights,
                                       std::string &consensus_sequence,
                                       std::string &consensus_sequence_gapped);


        /**
         * @brief Print in graphviz format.
         *
         */
        void print_as_dot(std::ostream &out_stream);


    private:
        std::unique_ptr<TheseusAlignerImpl> msa_aligner_impl_;
        bool still_end_to_end = true;
    };

} // namespace theseus
