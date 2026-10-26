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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stack>
#include <functional>
#include <utility>

#include "theseus/graph.h"
#include "theseus/alignment.h"

// Specific POA graph classes. Only used for MSA.
namespace theseus {

    using NodeId       = Graph::NodeId;
    using NodeView     = Graph::NodeView;
    using NodeIdRange  = Graph::NodeIdRange;
    using SequenceView = Graph::SequenceView;

    class POAVertex {
        public:
            std::vector<int> sequence_IDs;      // Sequence IDs
            std::vector<int> associated_vtxs;   // Associated vertexes
            std::vector<int> in_edges;          // In-going vertices
            std::vector<int> out_edges;         // Out-going vertices
            NodeId associated_node_compact;     // Corresponding node in the compact G graph
            int weight;                         // Weight of the node
            char value;                         // Base pair in this vertex
    };

    class POAEdge {
        public:
            int source;                         // Source vertex
            int destination;                    // Destination vertex
    };

    class POAGraph {
        public:

        std::vector<POAVertex> _poa_vertices;
        std::vector<POAEdge> _poa_edges;
        std::vector<int> _first_poa_vtx; // For each vertex in the compacted graph, the first POA vertex associated to it
        std::vector<int> _seq_weights;
        std::vector<int> _seq_starts;
        std::vector<int> _seq_ends;
        NodeId _end_vtx_poa;


        /**
         * @brief Split vertices in compacted_G if an edge splitting them is added.
         *
         * @param poa_source
         * @param poa_destination
         * @param compacted_G
         */
        void split_vertices(
            int poa_source,
            int poa_destination,
            Graph &compacted_G)
        {
            // Detect source and destination vertices in the compacted graph
            NodeId source_v_compact = _poa_vertices[poa_source].associated_node_compact;
            NodeId destination_v_compact = _poa_vertices[poa_destination].associated_node_compact;
            // Detect positions and necessity of split
            int pos_source = poa_source - _first_poa_vtx[source_v_compact];
            bool split_source = pos_source < compacted_G.node_size(source_v_compact) - 1;
            int pos_destination = poa_destination - _first_poa_vtx[destination_v_compact];
            bool split_destination = pos_destination > 0;
            if (source_v_compact == destination_v_compact && poa_source + 1 == poa_destination) {
                return; // No need to split
            }
            // Update source node
            if (split_source) {
                // New node
                int original_size  = compacted_G.node_size(source_v_compact);
                std::string suffix = compacted_G.split_sequence(source_v_compact, pos_source + 1);
                NodeId new_node_id = compacted_G.add_node(std::move(suffix));
                // Update the new node
                compacted_G.substitute_out_edges(source_v_compact, new_node_id);
                _first_poa_vtx.push_back(poa_source + 1);
                // Update the POA graph
                for (int l = pos_source + 1; l < original_size; ++l) {
                    _poa_vertices[(poa_source - pos_source) + l].associated_node_compact = new_node_id;
                }
                // Add the new edge
                compacted_G.add_edge(source_v_compact, new_node_id);
                // In case of self-node splitting gap, update destination_vtx
                if (source_v_compact == destination_v_compact) {
                    destination_v_compact = new_node_id;
                    pos_destination = poa_destination - _first_poa_vtx[new_node_id];
                }
            }
            // Update destination vertex
            if (split_destination) {
                // New node
                int original_size  = compacted_G.node_size(destination_v_compact);
                std::string suffix = compacted_G.split_sequence(destination_v_compact, pos_destination);
                NodeId new_node_id = compacted_G.add_node(std::move(suffix));
                // Update the new node
                compacted_G.substitute_out_edges(destination_v_compact, new_node_id);
                _first_poa_vtx.push_back(poa_destination);
                // Update the POA graph
                for (int l = pos_destination; l < original_size; ++l) {
                    _poa_vertices[(poa_destination - pos_destination) + l].associated_node_compact = new_node_id;
                }
                // Add the new edge
                compacted_G.add_edge(destination_v_compact, new_node_id);
                // New destination vertex
                destination_v_compact = new_node_id;
            }
            // Add the new edge
            compacted_G.add_edge(source_v_compact, destination_v_compact);
        }


        /**
         * @brief Update a vertex in the POA graph.
         *
         * @param poa_v             Vertex to update
         * @param new_node_id       Position of the new node in the compacted graph
         * @param value             Value to add
         * @param new_node_exists   Whether the node already exists or not in compacted_G
         * @param compacted_G       The compacted graph
         * @param seq_ID
         * @param added_weight
         */
        void update_poa_vertex(
            int &poa_v,
            NodeId &new_node_id,
            char value,
            bool &new_node_exists,
            Graph &compacted_G,
            int seq_ID,
            int added_weight)
        {
            // Check if the vertex already exists
            bool already_exists = false;
            int vtx;
            char vtx_value;
            for (long unsigned int l = 0; l < _poa_vertices[poa_v].associated_vtxs.size(); ++l)
            {
                vtx       = _poa_vertices[poa_v].associated_vtxs[l];
                vtx_value = _poa_vertices[vtx].value;
                if (value == vtx_value)
                {
                    already_exists    = true;
                    new_node_exists   = false;
                    poa_v = vtx; // poa_v is the vertex that will be used when adding an edge
                }
            }

            // Create it if it doesn't
            if (!already_exists)
            {
                POAVertex new_vertex;
                new_vertex.value           = value;
                new_vertex.associated_vtxs = _poa_vertices[poa_v].associated_vtxs;      // Associate it the necessary vertices
                new_vertex.associated_vtxs.push_back(poa_v);                            // Add the vtx poa_v, as it is missing
                new_vertex.weight          = 0;
                _poa_vertices.push_back(new_vertex);
                poa_v = _poa_vertices.size() - 1; // poa_v is the vertex that will be used when adding an edge

                // Update the other vertices
                for (long unsigned int l = 0; l < new_vertex.associated_vtxs.size(); ++l)
                {
                    vtx = new_vertex.associated_vtxs[l];
                    _poa_vertices[vtx].associated_vtxs.push_back(poa_v);
                }

                // Add a new vertex to the compacted graph
                if (new_node_exists) {
                    // Add a character to the existing vertex (the last one in compacted_G)
                    compacted_G.expand_sequence(new_node_id, std::string(1, value));
                    _poa_vertices[poa_v].associated_node_compact = new_node_id;
                }
                else {
                    // Create a new vertex (the last one in compacted_G)
                    new_node_id = compacted_G.add_node(std::string(1, value));
                    _first_poa_vtx.push_back(poa_v);
                    // Update the poa graph
                    _poa_vertices[poa_v].associated_node_compact = new_node_id;
                    new_node_exists = true;
                }
            }

            // Add sequence_ID
            _poa_vertices[poa_v].sequence_IDs.push_back(seq_ID);
            _poa_vertices[poa_v].weight += added_weight;
        }


        /**
         * @brief Update a POA edge.
         *
         * @param source
         * @param destination
         * @param G
         */
        void update_poa_edge(
            int source,
            int destination,
            Graph &compacted_G)
        {
            // Stop if one of the nodes is not valid
            if (source == -1 || destination == -1) return;

            // Check if the edge already exists
            bool already_exists = false;
            int curr_edge;
            for (long unsigned int l = 0; l < _poa_vertices[source].out_edges.size(); ++l) {
                curr_edge = _poa_vertices[source].out_edges[l];
                if (_poa_edges[curr_edge].source == source && _poa_edges[curr_edge].destination == destination) {
                    already_exists = true;
                }
            }
            // It doesn't, so you should create it
            if (!already_exists) {
                POAEdge new_edge;
                new_edge.source = source;
                new_edge.destination = destination;
                // Update the data in the poa graph
                _poa_vertices[source].out_edges.push_back(_poa_edges.size());
                _poa_vertices[destination].in_edges.push_back(_poa_edges.size());
                _poa_edges.push_back(new_edge);
                // Update the data in the compacted graph
                split_vertices(source, destination, compacted_G);
            }
        }

        void convert_path(
            Alignment &backtrace,
            std::vector<int> &poa_path,
            Graph &compacted_G,
            int end_column
        ) {
            // First vertex (add it because it is empty)
            poa_path.push_back(_first_poa_vtx[backtrace.path[0]]);
            // Convert the path to a path in the poa graph
            for (long unsigned int l = 1; l < backtrace.path.size() - 1; ++l) {
                int first_poa_vtx = _first_poa_vtx[backtrace.path[l]];
                // Add all internal vertices
                int node_size = compacted_G.node_size(backtrace.path[l]);
                for (int k = 0; k < node_size; ++k) {
                    poa_path.push_back(first_poa_vtx + k);
                }
            }
            // Add the nodes of last vertex (if it is not empty) and only up to the
            // end column.
            int last_vertex = backtrace.path[backtrace.path.size() - 1];
            int first_poa_vtx = _first_poa_vtx[last_vertex];
            for (int k = 0; k < end_column; ++k) {
                poa_path.push_back(first_poa_vtx + k);
            }
            // Connect to sink
            poa_path.push_back(_end_vtx_poa);
        }


        // Add alignment data in the _poa_graph
        void add_alignment_poa(
            Graph &compacted_G,
            Alignment &backtrace,
            SequenceView new_seq,
            int seq_ID,
            int weight,
            int end_column
        ) {
            // Convert the path to the corresponding path in the poa graph
            std::vector<int> poa_path;
            convert_path(backtrace, poa_path, compacted_G, end_column);
            bool new_node_exists = false;
            NodeId new_node_id;
            long unsigned int k = 0;
            int i = 0, l = 0, prev_v_poa = poa_path[0], new_v_poa = poa_path[0];
            while (k < backtrace.edit_op.size()) {
                if (backtrace.edit_op[k] == 'M') {  // Match
                    prev_v_poa = new_v_poa;
                    new_v_poa = poa_path[l + 1];
                    _poa_vertices[new_v_poa].sequence_IDs.push_back(seq_ID);
                    _poa_vertices[new_v_poa].weight += weight;
                    update_poa_edge(prev_v_poa, new_v_poa, compacted_G);
                    i += 1;
                    l += 1;
                    new_node_exists = false;
                }
                else if (backtrace.edit_op[k] == 'X') { // Mismatch
                    prev_v_poa = new_v_poa;
                    new_v_poa = poa_path[l + 1];
                    update_poa_vertex(new_v_poa, new_node_id, new_seq[i], new_node_exists, compacted_G, seq_ID, weight);
                    update_poa_edge(prev_v_poa, new_v_poa, compacted_G);
                    i += 1;
                    l += 1;
                }
                else if (backtrace.edit_op[k] == 'D') { // Deletion
                    // Add the new vertex
                    POAVertex new_vertex;
                    new_vertex.value  = new_seq[i];
                    new_vertex.sequence_IDs.push_back(seq_ID);
                    new_vertex.weight = weight;
                    _poa_vertices.push_back(new_vertex);
                    // Add/update a new compacted node
                    if (new_node_exists) {
                        // Add a character to the existing vertex (the last one in compacted_G)
                        compacted_G.expand_sequence(new_node_id, std::string(1, new_seq[i]));
                        _poa_vertices[_poa_vertices.size() - 1].associated_node_compact = new_node_id;
                    }
                    else {
                        // Add node
                        new_node_id = compacted_G.add_node(std::string(1, new_seq[i]));
                        _first_poa_vtx.push_back(_poa_vertices.size() - 1);
                        // Update the poa graph
                        _poa_vertices[_poa_vertices.size() - 1].associated_node_compact = new_node_id;
                        new_node_exists = true;
                    }
                    // Add the new edge
                    prev_v_poa = new_v_poa;
                    new_v_poa  = _poa_vertices.size() - 1;
                    update_poa_edge(prev_v_poa, new_v_poa, compacted_G);
                    i += 1;
                }
                else {
                    l += 1;
                }
                k += 1;
            }
            prev_v_poa = new_v_poa;
            // Add edge to the sink/source node, as no edit operation covers it?
            update_poa_edge(prev_v_poa, poa_path[poa_path.size() - 1], compacted_G); // Add the edge to the sink node
            // Update sequence weight, start poa vertex and end poa vertex
            _seq_weights.push_back(weight);
            if (new_seq.size() > 0) {
                _seq_starts.push_back(poa_path[1]);
                _seq_ends.push_back(prev_v_poa);
            }
            else {
                _seq_starts.push_back(-1);
                _seq_ends.push_back(-1);
            }
            // print_poa_dot();
        }


        /**
         * @brief Create the initial POA graph from the initial graph G.
         *
         */
        void create_initial_graph(theseus::Graph &G, int initial_weight)
        {
            // Source vertex
            theseus::POAVertex source_v;
            _first_poa_vtx.push_back(0);
            source_v.out_edges.push_back(0);
            source_v.associated_node_compact = 0;
            source_v.weight = initial_weight;
            source_v.value = '-';
            _poa_vertices.push_back(source_v);
            theseus::POAEdge source_edge;
            source_edge.source = 0;
            source_edge.destination = 1;
            _poa_edges.push_back(source_edge);
            // Central vertices
            NodeView node_view = G.node(1);
            _first_poa_vtx.push_back(1);
            for (int l = 0; l < G.node_size(1); ++l) {
                theseus::POAVertex new_v;
                new_v.in_edges.push_back(_poa_edges.size() - 1);
                new_v.out_edges.push_back(_poa_edges.size());
                new_v.value = node_view.sequence[l];
                new_v.associated_node_compact = 1;
                new_v.weight = initial_weight;
                new_v.sequence_IDs.push_back(0); // Sequence ID 0
                _poa_vertices.push_back(new_v);
                theseus::POAEdge new_edge;
                new_edge.source = _poa_vertices.size() - 1;
                new_edge.destination = _poa_vertices.size();
                _poa_edges.push_back(new_edge);
            }
            // Sink vertex
            theseus::POAVertex sink_v;
            _first_poa_vtx.push_back(G.node_size(1) + 1);
            sink_v.in_edges.push_back(_poa_edges.size() - 1);
            sink_v.associated_node_compact = 2;
            sink_v.weight = initial_weight;
            sink_v.value = '-';
            _poa_vertices.push_back(sink_v);
            // Set the end vertex of the POA graph
            _end_vtx_poa = _poa_vertices.size() - 1;

            // Update start poa vertex, end poa vertex and sequence weight
            if (G.node_size(1) > 0) {
                _seq_weights.push_back(initial_weight);
                _seq_starts.push_back(1);
                _seq_ends.push_back(G.node_size(1));
            }
            else {
                _seq_weights.push_back(initial_weight);
                _seq_starts.push_back(-1);
                _seq_ends.push_back(-1);
            }
        }


        /**
         * @brief Compute the consensus sequence and the MSA matrix of the POA
         * graph based on the majority voting algorithm.
         *
         * The majority voting algorithm works as follows:
         * 1) You topologically order the vertices of the augmented POA graph (with
         * extra edges between aligned nodes) using DFS.
         * 2) You assign a column index to each vertex in the MSA representation.
         * 3) During alignment, you have stored the start and end poa vertices
         * for each sequence, along with its weight.
         * 4) With these values, you can determine, for each column, the wighted number
         * of sequences that are valid in that column. That is, the number of weighted
         * sequences whose interval of valid columns includes the column you are
         * looking at.
         * 4) The consensus sequence has a base pair per column where more than half
         * of the weighted sequences do not have a gap there. The base pair is
         * the one with the highest weight among the vertices in that column.
         *
         */
        void poa_fasta_and_majority(int num_sequences,
                                    std::string &consensus_sequence,
                                    std::string &consensus_sequence_gapped,
                                    std::vector<int> &consensus_weights,
                                    std::vector<std::vector<char>> &msa,
                                    size_t &source_id,
                                    size_t &sink_id) {

            // Create an augmented graph to ensure MSA integrity (to ensure valid topological order)
            POAGraph augmented_poa_graph;
            augmented_poa_graph._poa_vertices = _poa_vertices;
            augmented_poa_graph._poa_edges = _poa_edges;

            // Given an edge e = (v1, v2), add an extra edge per aliged pair of
            // the aligned nodes to v1 and v2. That is, new_e = (w1, w2) where
            // w1 and w2 are aligned nodes to v1 and v2, respectively.
            int num_original_edges = augmented_poa_graph._poa_edges.size();
            for (int l = 0; l < num_original_edges; ++l) {
                POAEdge edge = augmented_poa_graph._poa_edges[l];
                // For each edge, add an extra edge for each aligned pair
                POAVertex &source_vertex      = augmented_poa_graph._poa_vertices[edge.source];
                POAVertex &destination_vertex = augmented_poa_graph._poa_vertices[edge.destination];
                for (long unsigned int i = 0; i < source_vertex.associated_vtxs.size(); ++i) {
                    int aligned_source = source_vertex.associated_vtxs[i];
                    for (long unsigned int j = 0; j < destination_vertex.associated_vtxs.size(); ++j) {
                        int aligned_destination = destination_vertex.associated_vtxs[j];
                        // Create a new edge between the aligned nodes
                        POAEdge new_edge;
                        new_edge.source = aligned_source;
                        new_edge.destination = aligned_destination;
                        augmented_poa_graph._poa_edges.push_back(new_edge);
                        // Add the information to the vertices
                        augmented_poa_graph._poa_vertices[aligned_source].out_edges.push_back(augmented_poa_graph._poa_edges.size() - 1);
                        augmented_poa_graph._poa_vertices[aligned_destination].in_edges.push_back(augmented_poa_graph._poa_edges.size() - 1);
                    }
                }
                // Add source to aligned destination edges
                for (size_t j = 0; j < destination_vertex.associated_vtxs.size(); ++j) {
                    int aligned_destination = destination_vertex.associated_vtxs[j];
                    // Create a new edge between the source and the aligned destination node
                    POAEdge new_edge;
                    new_edge.source = edge.source;
                    new_edge.destination = aligned_destination;
                    augmented_poa_graph._poa_edges.push_back(new_edge);
                    // Add the information to the vertices
                    source_vertex.out_edges.push_back(augmented_poa_graph._poa_edges.size() - 1);
                    augmented_poa_graph._poa_vertices[aligned_destination].in_edges.push_back(augmented_poa_graph._poa_edges.size() - 1);
                }
                // Add aligned source to destination edges
                for (size_t j = 0; j < source_vertex.associated_vtxs.size(); ++j) {
                    int aligned_source = source_vertex.associated_vtxs[j];
                    // Create a new edge between the aligned source and the destination node
                    POAEdge new_edge;
                    new_edge.source = aligned_source;
                    new_edge.destination = edge.destination;
                    augmented_poa_graph._poa_edges.push_back(new_edge);
                    // Add the information to the vertices
                    augmented_poa_graph._poa_vertices[aligned_source].out_edges.push_back(augmented_poa_graph._poa_edges.size() - 1);
                    destination_vertex.in_edges.push_back(augmented_poa_graph._poa_edges.size() - 1);
                }
            }

            // Topologically order the vertices using DFS
            std::vector<bool> visited(augmented_poa_graph._poa_vertices.size(), false);
            std::stack<int> topo_stack;

            // Stacked dfs to perform the topological sort. It is important to
            // avoid a recursive dfs, as it causes stack overflow for large graphs.
            // Perform DFS for all vertices to handle disconnected components
            for (size_t start_v = 0; start_v < augmented_poa_graph._poa_vertices.size(); ++start_v) {
                if (!visited[start_v]) {
                    // Explicit stack instead of recursion. The stack stores
                    // pairs: {vertex_id, neighbors_processed_flag}
                    std::stack<std::pair<int, bool>> st;
                    st.push({start_v, false});

                    // Traversal from start_v
                    while (!st.empty()) {
                        auto [v, neighbors_processed] = st.top();
                        st.pop();

                        if (!neighbors_processed) {
                            // If it was already visited by another path while waiting in the stack, skip
                            if (visited[v]) {
                                continue;
                            }

                            // Mark as visited now that we are actively processing it
                            visited[v] = true;

                            // Push it back onto the stack with the flag set to 'true'.
                            // When we pop it again, we will know that all of its neighbors,
                            // located higher on the stack, have been processed and popped.
                            st.push({v, true});

                            // Push all unvisited neighbors onto the intermediate stack "st"
                            for (int edge_idx : augmented_poa_graph._poa_vertices[v].out_edges) {
                                int next_v = augmented_poa_graph._poa_edges[edge_idx].destination;
                                if (!visited[next_v]) {
                                    st.push({next_v, false});
                                }
                            }
                        } else {
                            // Popping a node whose nodes have been processed means
                            // that we can add it to the topological order stack.
                            topo_stack.push(v);
                        }
                    }
                }
            }

            // Reverse the stack to get the topological order
            std::vector<int> topological_order;
            while (!topo_stack.empty()) {
                topological_order.push_back(topo_stack.top());
                topo_stack.pop();
            }

            // Determine the columns of the nodes in the MSA representation
            std::vector<int> node_to_column;
            node_to_column.resize(augmented_poa_graph._poa_vertices.size(), -1);
            int column_index = 0;
            for (int v : topological_order) {
                // Check aligned nodes
                POAVertex &vertex = augmented_poa_graph._poa_vertices[v];
                int min_aligned = -1;
                for (int aligned_v : vertex.associated_vtxs) {
                    if (node_to_column[aligned_v] != -1) {
                        min_aligned = node_to_column[aligned_v];
                    }
                }

                if (min_aligned != -1) {
                    node_to_column[v] = min_aligned;  // Assign the column of the aligned node
                } else {
                    node_to_column[v] = column_index; // Assign a new column
                    column_index += 1;
                }
            }
            // Determine column to node mapping
            std::vector<std::vector<int>> column_to_nodes;
            column_to_nodes.resize(column_index);
            for (long unsigned int v = 0; v < node_to_column.size(); ++v) {
                int column = node_to_column[v];
                if (column != -1) {
                    column_to_nodes[column].push_back(v);
                }
            }

            // Create MSA matrix
            int columns = column_index;   // Number of columns in the MSA
            int rows = num_sequences + 2; // Number of sequences  (+1 for original node +1 for consensus sequence)
            msa.resize(rows, std::vector<char>(columns, '-'));

            // Fill the MSA with the aligned sequences (except first and last nodes)
            for (size_t l = 1; l < _poa_vertices.size(); ++l) {
                POAVertex &vertex = augmented_poa_graph._poa_vertices[l];
                int column = node_to_column[l];
                // Find the sequence IDs of the node
                for (size_t k = 0; k < vertex.sequence_IDs.size(); ++k) {
                    int seq_id = vertex.sequence_IDs[k];
                    // Fill the MSA with the value of the vertex
                    msa[seq_id][column] = vertex.value;
                }
            }
            source_id = node_to_column[0];
            sink_id   = node_to_column[_end_vtx_poa];

            ///////////////////// Majority voting consensus ////////////////////
            // Create two vectors indicating the weighted number of sequences that
            // start and end in each column, respectively
            std::vector<int> seq_starts_in_column(column_index, 0);
            std::vector<int> seq_ends_in_column(column_index, 0);
            // The data on the start and end poa vertices is in:
            //    - std::vector<int> _seq_starts;
            //    - std::vector<int> _seq_ends;
            // "num_sequences" is the ID of the last sequence, so there are
            // _seq_starts.size() = num_sequences + 1 sequences. Empty sequences
            // have no start/end vertices (-1).
            for (size_t i = 0; i < _seq_starts.size(); ++i) {
                if (_seq_starts[i] < 0) continue;
                int start_column = node_to_column[_seq_starts[i]];
                int end_column   = node_to_column[_seq_ends[i]];
                int weight       = _seq_weights[i];
                if (start_column != -1) {
                    seq_starts_in_column[start_column] += weight;
                }
                if (end_column != -1) {
                    seq_ends_in_column[end_column] += weight;
                }
            }
            // Compute the weighted number of sequences that are valid in each column
            std::vector<int> valid_sequences_in_column(column_index, 0);
            valid_sequences_in_column[0] = seq_starts_in_column[0];
            for (int j = 1; j < column_index; ++j) {
                valid_sequences_in_column[j] = valid_sequences_in_column[j - 1]
                + seq_starts_in_column[j]       // Newly activated
                - seq_ends_in_column[j - 1];    // No longer valid
            }

            // Construct the consensus sequence based on the weighted majority voting
            consensus_sequence = "";
            for (int j = 0; j < column_index; ++j) {
                // Count the number of sequences that do not have a gap in that column
                int non_gap_weight = 0;
                for (int v : column_to_nodes[j]) {
                    non_gap_weight += augmented_poa_graph._poa_vertices[v].weight;
                }
                // Check if it is more than half of the valid sequences in that column
                // std::cout << "Column " << j << ": " << non_gap_weight << " non-gap weight, " << valid_sequences_in_column[j] << " valid sequences\n";
                if (non_gap_weight > valid_sequences_in_column[j] / 2) {
                    int curr_max_weight = -1;
                    char curr_consensus_char = '-';
                    // Find the bp with highest weight
                    for (int v : column_to_nodes[j]) {
                        // std::cout << "Vertex " << v << ": value = " << augmented_poa_graph._poa_vertices[v].value << ", weight = " << augmented_poa_graph._poa_vertices[v].weight << "\n";
                        if (augmented_poa_graph._poa_vertices[v].weight > curr_max_weight) {
                            curr_max_weight = augmented_poa_graph._poa_vertices[v].weight;
                            curr_consensus_char = augmented_poa_graph._poa_vertices[v].value;
                        }
                    }
                    if (curr_consensus_char != '-') {
                        consensus_sequence += curr_consensus_char;
                        consensus_weights.push_back(curr_max_weight);
                    }
                    // Add the consensus sequence character to the MSA matrix
                    msa[msa.size() - 1][j] = curr_consensus_char;
                }
            }

            // The last row of the MSA matrix is the consensus sequence with gaps
            // Convert from vector<char> to string
            for (size_t j = 0; j < msa[msa.size() - 1].size(); ++j) {
                if (j != source_id && j != sink_id) {
                    consensus_sequence_gapped += msa[msa.size() - 1][j];
                }
            }
        }


        /**
         * @brief Convert the POA graph to a FASTA file (MSA format)
         *
         * @param output_file
         */
        void poa_to_fasta(int num_sequences, std::ostream &out_file) {
            // Get column ordering and nodes in each column
            std::vector<int> consensus_weights;
            std::vector<std::vector<char>> msa;
            std::string consensus_sequence;
            std::string consensus_sequence_gapped;
            size_t source_id, sink_id;
            poa_fasta_and_majority(num_sequences, consensus_sequence, consensus_sequence_gapped, consensus_weights, msa, source_id, sink_id);

            // Print the MSA
            for (size_t i = 0; i < msa.size(); ++i) {
                if (i < msa.size() - 1) {
                    out_file << ">Sequence_" << i + 1 << "\n"; // Sequence ID
                } else {
                    out_file << ">Consensus\n";
                }
                for (size_t j = 0; j < msa[i].size(); ++j) {
                    if (j != source_id && j != sink_id) {
                        out_file << msa[i][j];
                    }
                }
                out_file << "\n"; // New line after each sequence
            }
        }


        /**
         * @brief Print the consensus sequence of the POA graph
         *
         */
        std::string poa_to_consensus() {
            // Apply heaviest bundle algorithm to find the consensus sequence
            std::string consensus_sequence = "";
            NodeId current_vertex = 0; // Start from the source vertex

            while (current_vertex != _end_vtx_poa) {
                POAVertex &vertex = _poa_vertices[current_vertex];

                // Find the outgoing edge with the maximum weight
                int max_weight = -1;
                int next_vertex = -1;
                for (int edge_idx : vertex.out_edges) {
                    POAEdge &edge = _poa_edges[edge_idx];
                    int destination_vtx = edge.destination;
                    if (_poa_vertices[destination_vtx].weight > max_weight) {
                        max_weight = _poa_vertices[destination_vtx].weight;
                        next_vertex = destination_vtx;
                    }
                }

                // If no outgoing edges, break the loop (should't happen in a well-formed POA graph)
                if (next_vertex == -1) {
                    std::cerr << "No outgoing edges from vertex " << current_vertex << ". Ending consensus extraction.\n";
                    break;
                }

                // Append the value of the next vertex to the consensus sequence
                consensus_sequence += _poa_vertices[next_vertex].value;

                // Move to the next vertex
                current_vertex = next_vertex;
            }

            // Elimimate last character, corresponding to the sink node
            if (!consensus_sequence.empty()) {
                consensus_sequence.pop_back();
            }
            return consensus_sequence;
        }

        /**
         * @brief return the consensus sequence of the POA graph based on the majority voting algorithm
         *
         */
        void poa_to_consensus_weighted_majority_voting(int num_sequences,
                                                       std::vector<int> &consensus_weights,
                                                       std::string &consensus_sequence,
                                                       std::string &consensus_sequence_gapped) {
            size_t source_id, sink_id;
            std::vector<std::vector<char>> msa_matrix;
            // Clear previous data if there was some
            consensus_sequence = "";
            consensus_sequence_gapped = "";
            consensus_weights.clear();
            poa_fasta_and_majority(num_sequences, consensus_sequence, consensus_sequence_gapped, consensus_weights, msa_matrix, source_id, sink_id);
        }


        /*
         * Print the POA graph in DOT format (for debugging purposes)
        */
        void print_poa_dot() {
            std::cout << "digraph POA {\n";
            for (size_t i = 0; i < _poa_vertices.size(); ++i) {
                // Print weight
                std::cout << "  " << i << " [label=\"" << _poa_vertices[i].value << " " << _poa_vertices[i].weight << "\"];\n";
            }
            for (size_t i = 0; i < _poa_edges.size(); ++i) {
                int source = _poa_edges[i].source;
                int destination = _poa_edges[i].destination;
                std::cout << "  " << source << " -> " << destination << ";\n";
            }
            std::cout << "}\n";
        }
    };
}