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
#include <vector>
#include "theseus/graph.h"
#include "theseus/penalties.h"

/**
 * Class containing the alignment information (CIGAR and path) and a custom conversion
 * score according to user defined penalties.
 *
 */

namespace theseus {

using NodeId = Graph::NodeId;

/**
 * Alignment modes:
 * - Full: compute the score and the alignment (CIGAR and path). The M
 *   wavefronts are kept until the end of the alignment for the backtrace.
 * - ScoreOnly: compute only the score and the end position of the alignment.
 *   Only the wavefronts of the last scores (the scope) are kept in memory.
 */
enum class AlignmentMode {
    Full,
    ScoreOnly
};

/**
 * End conditions of an alignment to the POA graph of the MSA. The aligner stops
 * as soon as the condition is met:
 * - FullQuery (default): the whole query is consumed. The alignment can end at
 *   any position of the graph. The sequence is always linked to the sink (or to
 *   the source in a reverse alignment) in the POA graph.
 * - EndToEnd: the whole query is consumed AND the end of the graph (sink, or
 *   source in a reverse alignment) is reached.
 * - EndsFree: the whole query is consumed OR the end of the graph is reached.
 *   If the end of the graph is not reached, the sequence is not linked to it.
 */
enum class EndCondition {
    FullQuery,
    EndsFree,
    EndToEnd
};

/**
 * Options of the alignment of one sequence to the POA graph of the MSA:
 * - lag_pruning:   prune the diagonals that lag behind in the alignment.
 * - reverse:       align from the end of the graph (the sink) backwards.
 * - end_condition: see EndCondition.
 * - density_drop:  drop the alignment when it advances too slowly. The partial
 *                  alignment found is added to the POA graph and the returned
 *                  status is THESEUS_STATUS_END_UNREACHABLE.
 */
struct MSAOptions {
    bool         lag_pruning   = false;
    bool         reverse       = false;
    EndCondition end_condition = EndCondition::FullQuery;
    bool         density_drop  = false;
};

class Alignment {
    public:
    /*
    * Error codes & messages
    */
    // [OK]
    #define THESEUS_STATUS_ALG_COMPLETED            0  // Success (Complete alignment found)
    #define THESEUS_STATUS_ALG_PARTIAL              1  // Success (Partial alignment found)
    // [FAIL]
    #define THESEUS_STATUS_MAX_STEPS_REACHED     -100  // Maximum number of Theseus-steps reached
    #define THESEUS_STATUS_UNATTAINABLE          -300  // Alignment unattainable under configured heuristics
    // [INTERNAL]
    #define THESEUS_STATUS_OK                      -1  // Computing alignment (in progress)
    #define THESEUS_STATUS_END_REACHED             -2  // Alignment end reached
    #define THESEUS_STATUS_END_UNREACHABLE         -3  // Alignment end unreachable under current configuration (eg advancement density)


    /*
    * Error messages
    */
    // OK
    #define THESEUS_STATUS_ALG_COMPLETED_MSG           "[Theseus] Alignment completed successfully"
    #define THESEUS_STATUS_ALG_PARTIAL_MSG             "[Theseus] Alignment extension computed (partial alignment)"
    #define THESEUS_STATUS_ALG_COMPLETED_MSG_SHORT     "OK.Full"
    #define THESEUS_STATUS_ALG_PARTIAL_MSG_SHORT       "OK.Partial"
    // FAILED
    #define THESEUS_STATUS_MAX_STEPS_REACHED_MSG       "[Theseus] Alignment failed. Maximum Theseus-steps limit reached"
    #define THESEUS_STATUS_UNATTAINABLE_MSG            "[Theseus] Alignment failed. Unattainable under configured heuristics"
    #define THESEUS_STATUS_MAX_STEPS_REACHED_MSG_SHORT "FAILED.MaxTheseusSteps"
    #define THESEUS_STATUS_UNATTAINABLE_MSG_SHORT      "FAILED.Unattainable"
    // Internal
    #define THESEUS_STATUS_END_REACHED_MSG             "[Theseus] Alignment end reached"
    #define THESEUS_STATUS_END_UNREACHABLE_MSG         "[Theseus] Alignment end unreachable under current configuration (due to heuristics)"
    #define THESEUS_STATUS_UNKNOWN_MSG                 "[Theseus] Unknown error code"
    #define THESEUS_STATUS_END_REACHED_MSG_SHORT       "INTERNAL.Reached"
    #define THESEUS_STATUS_END_UNREACHABLE_MSG_SHORT   "INTERNAL.Dropped"
    #define THESEUS_STATUS_UNKNOWN_MSG_SHORT           "Unknown"


      /**
       * @brief Backtrace objects: similar to the CIGAR in sequence alignment. It consists
       *        of the set of edit operations of the alignment and the path of the alignment
       *        through the reference graph.
       */
      std::vector<char> edit_op;    // Edit operations (empty in score-only mode)
      std::vector<NodeId> path;     // Path of the alignment (in score-only mode,
                                    // only the last vertex of the path)
      int start_offset;             // Start offset in the first vertex of the path
      int end_offset;               // End offset in the last vertex of the path
      int theseus_status;           // Alignment status
      int score = 0;                // Alignment score according to the user penalties.
                                    // Partial alignments: score of the partial path
                                    // (consumed prefix of the sequence)

    private:
};

} // namespace theseus