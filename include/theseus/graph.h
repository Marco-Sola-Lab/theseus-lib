#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <iostream>

/**
 * The internal graph used by the Theseus aligner.
 *
 * The graph supports adding and removing nodes and edges. Each node of the graph
 * is identified by an id that is stable until the node is removed or the graph
 * is destroyed. The graph also supports retrieving a non-owning view of a node,
 * which includes the node's sequence and its incoming and outgoing edges.
 *
 * WARNING: Any modification to the graph (adding or removing nodes or edges)
 * may invalidate existing views and iterators.
 *
 */

namespace theseus {

class GraphImpl;

class Graph {
public:
    // Numeric identifier for a node.
    using NodeId = std::size_t;

    /**
     * Exception thrown when a NodeId is requested that does not exist
     * within the graph's current state.
     */
    class InvalidNodeException : public std::out_of_range {
    public:
        explicit InvalidNodeException(std::size_t id)
            : std::out_of_range("NodeId " + std::to_string(id) +
                                " is not present in the graph."),
              id_(id) {}

        [[nodiscard]] std::size_t id() const noexcept { return id_; }

    private:
        std::size_t id_;
    };

    /**
     * Range of NodeIds compliant with C++20std::ranges::range concept.
     *
     * WARNING: Ranges and internal iterators are invalidated when the graph
     * is modified.
     */
    class NodeIdRange {
    public:
        /**
         * Forward iterator for NodeIdRange. The concrete iterator type is erased
         * and stored in an inline buffer to keep implementation details hidden
         * and avoid dynamic allocations for performance.
         *
         */
        class iterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using value_type = NodeId;
            using difference_type = std::ptrdiff_t;
            using pointer = const NodeId *;
            using reference = const NodeId &;

            iterator() = default;

            ~iterator() { destroy(); }

            iterator(const iterator &o) { copy_from(o); }

            iterator &operator=(const iterator &o) {
                if (this != &o) {
                    destroy();
                    copy_from(o);
                }
                return *this;
            }

            iterator(iterator &&o) noexcept { move_from(o); }

            iterator &operator=(iterator &&o) noexcept {
                if (this != &o) {
                    destroy();
                    move_from(o);
                }
                return *this;
            }

            reference operator*() const { return *ops_->deref(buff_); }

            pointer operator->() const { return ops_->deref(buff_); }

            iterator &operator++() {
                ops_->incr(buff_);
                return *this;
            }

            iterator operator++(int) {
                auto c = *this;
                ++*this;
                return c;
            }

            bool operator==(const iterator &o) const noexcept {
                if (ops_ != o.ops_) {
                    return false;
                }
                if (ops_ == nullptr) {
                    return true;
                }
                return ops_->equal(buff_, o.buff_);
            }

            bool operator!=(const iterator &o) const noexcept { return !(*this == o); }

        private:
            // Only NodeIdRange can construct iterators.
            friend class NodeIdRange;

            // Size of the inline buffer for storing the concrete iterator.
            // This should be large enough to hold the largest concrete iterator
            // type.
            static constexpr std::size_t buff_size_bytes = 24;

            // Ops for type-erased iterator. This is a struct that contains
            // function pointers for the required iterator operations.
            struct Ops {
                const NodeId *(*deref)(const std::byte *);
                void (*incr)(std::byte *);
                bool (*equal)(const std::byte *, const std::byte *);
                void (*copy)(std::byte *dst, const std::byte *src);
                void (*move)(std::byte *dst, std::byte *src);
                void (*dtor)(std::byte *);
            };

            alignas(std::max_align_t) std::byte buff_[buff_size_bytes]{};
            const Ops *ops_ = nullptr;

            /**
             * Construct an iterator from a concrete iterator type. The concrete
             * iterator must satisfy the requirements of a forward iterator and
             * its size must not exceed buff_size_bytes. The concrete iterator
             * is type-erased and stored in the inline buffer, and the
             * corresponding operations are stored in the Ops struct.
             *
             * @tparam It The type of the concrete iterator.
             * @param it The concrete iterator.
             */
            template <typename It>
            explicit iterator(It it) {
                static_assert(
                    sizeof(It) <= buff_size_bytes,
                    "Concrete iterator exceeds inline buffer. Increase buff_size_bytes.");

                static constexpr Ops kops {
                    [](const std::byte *b) -> const NodeId * {
                        return std::addressof(*(*reinterpret_cast<const It *>(b)));
                    },
                    [](std::byte *b) {
                        ++(*reinterpret_cast<It *>(b));
                    },
                    [](const std::byte *a, const std::byte *b) -> bool {
                        return *reinterpret_cast<const It *>(a) ==
                               *reinterpret_cast<const It *>(b);
                    },
                    [](std::byte *dst, const std::byte *src) {
                        new (dst) It(*reinterpret_cast<const It *>(src));
                    },
                    [](std::byte *dst, std::byte *src) {
                        new (dst) It(std::move(*reinterpret_cast<It *>(src)));
                    },
                    [](std::byte *b) {
                        reinterpret_cast<It *>(b)->~It();
                    },
                };

                new (buff_) It(std::move(it));
                ops_ = &kops;
            }

            void destroy() {
                if (ops_) {
                    ops_->dtor(buff_);
                    ops_ = nullptr;
                }
            }

            void copy_from(const iterator &o) {
                ops_ = o.ops_;
                if (ops_)
                    ops_->copy(buff_, o.buff_);
            }

            void move_from(iterator &o) noexcept {
                ops_ = o.ops_;
                if (ops_) {
                    ops_->move(buff_, o.buff_);
                    o.ops_ = nullptr;
                }
            }
        };

        template <typename It>
        NodeIdRange(It begin, It end) : begin_(begin), end_(end) {}

        [[nodiscard]] iterator begin() const { return begin_; }
        [[nodiscard]] iterator end() const { return end_; }

    private:
        iterator begin_;
        iterator end_;
    };

    class SequenceView {
    public:
        struct iterator {
            using iterator_category = std::random_access_iterator_tag;
            using value_type = char;
            using difference_type = std::ptrdiff_t;
            using pointer = const char *;
            using reference = const char &;

            const char *ptr;
            int step;  // +1 forward, -1 reverse
            int deref_offset;  // 0 for forward, -1 for reverse

            reference operator[](difference_type n) const {
                return *(*this + n);
            }

            reference operator*() const {
                return *(ptr + deref_offset);
            }

            pointer operator->() const {
                return ptr + deref_offset;
            }

            iterator &operator++() {
                ptr += step;
                return *this;
            }

            iterator operator++(int) {
                iterator tmp = *this;
                ptr += step;
                return tmp;
            }

            iterator &operator--() {
                ptr -= step;
                return *this;
            }

            iterator operator--(int) {
                iterator tmp = *this;
                ptr -= step;
                return tmp;
            }

            iterator &operator+=(difference_type n) {
                ptr += n * step;
                return *this;
            }

            iterator &operator-=(difference_type n) {
                ptr -= n * step;
                return *this;
            }

            iterator operator+(difference_type n) const {
                return {ptr + n * step, step, deref_offset};
            }

            iterator operator-(difference_type n) const {
                return {ptr - n * step, step, deref_offset};
            }

            friend iterator operator+(difference_type n, const iterator &it) {
                return it + n;
            }

            difference_type operator-(const iterator &other) const {
                return (ptr - other.ptr) * step;
            }

            bool operator==(const iterator &o) const {
                return ptr == o.ptr && step == o.step;
            }

            bool operator!=(const iterator &o) const { return !(*this == o); }

            bool operator<(const iterator &o) const { return (ptr - o.ptr) * step < 0; }
            bool operator>(const iterator &o) const { return (ptr - o.ptr) * step > 0; }
            bool operator<=(const iterator &o) const { return (ptr - o.ptr) * step <= 0; }
            bool operator>=(const iterator &o) const { return (ptr - o.ptr) * step >= 0; }
        };

        SequenceView(std::string_view sv, bool reversed = false)
            : data_(sv), reversed_(reversed) {}

        iterator begin() const {
            if (!reversed_) {
                return {data_.data(), 1, 0};
            }
            else {
                return {data_.data() + data_.size(), -1, -1};
            }
        }

        iterator end() const {
            if (!reversed_) {
                return {data_.data() + data_.size(), 1, 0};
            }
            else {
                return {data_.data(), -1, -1};
            }
        }

        size_t size() const { return data_.size(); }

        // Underlying characters, always in forward order (data()[size() - 1 - idx]
        // is the character idx of a reversed view)
        const char *data() const { return data_.data(); }

        bool empty() const { return data_.empty(); }

        char operator[](size_t idx) const {
            return reversed_ ? data_[data_.size() - 1 - idx] : data_[idx];
        }

        // Output
        friend std::ostream& operator<<(std::ostream& os, const SequenceView& sv) {
            for (char c : sv) {
                os << c;
            }
            return os;
        }

        void change_reversed_flag(bool new_flag_value) {
            reversed_ = new_flag_value;
        }

        [[nodiscard]] bool is_reversed() const noexcept { return reversed_; }

        [[nodiscard]] SequenceView reversed() const noexcept {
            return SequenceView(data_, !reversed_);
        }

    private:
        std::string_view data_;
        bool reversed_;
    };


    /**
     * A non-owning view of a node in the graph. The view includes the node's
     * sequence and its incoming and outgoing edges.
     *
     * WARNING: Views are invalidated when the graph is modified.
     */
    struct NodeView {
    public:
        NodeView(SequenceView sequence_, NodeIdRange in_nodes_, NodeIdRange out_nodes_)
            : sequence(sequence_), in_nodes(in_nodes_), out_nodes(out_nodes_) {}

        SequenceView sequence;
        NodeIdRange in_nodes;
        NodeIdRange out_nodes;
    };

    /**
     * Construct an empty graph.
     *
     */
    Graph();

    /**
     * Destroy the graph and all its nodes and edges.
     *
     */
    ~Graph();

    /**
     * Copy constructor. Creates a deep copy of the graph.
     *
     * @param other The graph to copy.
     */
    Graph(const Graph& other);

    /**
     * Copy assignment operator. Creates a deep copy of the graph.
     *
     * @param other The graph to copy.
     * @return A reference to the current graph.
     */
    Graph& operator=(const Graph& other);

    /**
     * Move constructor.
     *
     * @param other The graph to move.
     */
    Graph(Graph&&) noexcept;

    /**
     * Move assignment operator.
     *
     * @param other The graph to move.
     * @return A reference to the current graph.
     */
    Graph& operator=(Graph&&) noexcept;

    /**
     * Swap the contents of this graph with another graph.
     *
     * @param other The graph to swap with.
     */
    void swap(Graph& other) noexcept;

    /**
     * Add a node with a copy of the given sequence to the graph and return its
     * id. The new node has no incoming or outgoing edges.
     *
     * @param sequence The sequence to associate with the new node.
     * @return The id of the newly added node.
     */
    [[nodiscard]] NodeId add_node(std::string_view sequence);

    /**
     * Add a node from a null-terminated C string.
     *
     * This overload exists to avoid ambiguity between std::string_view and
     * std::string&& overloads when passing string literals.
     *
     * @param sequence The null-terminated C string to associate with the new node.
     * @return The id of the newly added node.
     */
    [[nodiscard]] NodeId add_node(const char *sequence);

    /**
     * Add a node with the given sequence to the graph and return its id.
     * The new node has no incoming or outgoing edges.
     *
     * @param sequence The sequence to associate with the new node. The contents of
     * the string are moved into the graph.
     * @return The id of the newly added node.
     */
    [[nodiscard]] NodeId add_node(std::string&& sequence);

    /**
     * @brief Expand the sequence of a node with a suffix.
     *
     * @param id The id of the node.
     * @param suffix The suffix to append to the node's sequence.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    void expand_sequence(NodeId id, std::string_view suffix);

    /**
     * Split the sequence of a node at the given index, modifying the node
     * in-place and returning the detached suffix.
     *
     * The node retains the prefix [0, idx) and the returned string contains
     * the suffix [idx, len). Splitting at 0 leaves the node's sequence empty
     * and returns the full sequence; splitting at len leaves the node unchanged
     * and returns an empty string.
     *
     * @param id The id of the node.
     * @param idx Split position in [0, len], where len is the length of the
     * node's sequence.
     * @return The suffix of the sequence starting at idx.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    std::string split_sequence(NodeId id, size_t idx);

    /**
     * Remove the node with the given id from the graph, along with all its
     * incoming and outgoing edges. After a node is removed, its id may be
     * reused for new nodes added to the graph.
     *
     * @param id The id of the node to remove.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    void remove_node(NodeId id);

    /**
     * Add a directed edge from one node to another. Duplicate edges (multiple
     * edges from the same source to the same target) are allowed
     *
     * @param from The id of the node from which the edge originates.
     * @param to The id of the node to which the edge points.
     * @throws InvalidNodeException if either NodeId does not correspond to an
     * existing node.
     */
    void add_edge(NodeId from, NodeId to);

    /**
     * Remove a directed edge from one node to another. If multiple edges exist
     * from the same source to the same target, only one of them is removed.
     *
     * @param from The id of the node from which the edge originates.
     * @param to The id of the node to which the edge points.
     * @return true if an edge was found and removed, false otherwise.
     * @throws InvalidNodeException if either NodeId does not correspond to an
     * existing node.
     */
    bool remove_edge(NodeId from, NodeId to);

    /**
     * @brief Remove all outgoing edges from the given node.
     *
     * @param from The id of the node from which to remove all outgoing edges.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    void remove_out_edges(NodeId from);

    /**
     * @brief Remove all incoming edges from the given node.
     *
     * @param from The id of the node from which to remove all incoming edges.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    void remove_in_edges(NodeId from);

    /**
     * @brief Given an "original_from" node and a "new_from" node, substitute
     * all outgoing edges of the "orig_from" node to originate from the "new_from"
     * node. This is used for example when splitting a node: the new node
     * created by the split now inherits the outgoing edges of the original
     * node.
     *
     * @param orig_from The id of the original from node.
     * @param new_from The id of the new from node.
     * @throws InvalidNodeException if either NodeId does not correspond to an
     * existing node.
     */
    void substitute_out_edges(NodeId orig_from, NodeId new_from);

    /**
     * @brief Given an "original_to" node and a "new_to" node, substitute
     * all incoming edges of orig_to node to point to new_to node. This is
     * used for example when splitting a node: the new node created by the
     * split now inherits the incoming edges of the original node.
     *
     * @param orig_to The id of the original to node.
     * @param new_to The id of the new to node.
     * @throws InvalidNodeException if either NodeId does not correspond to an
     * existing node.
     */
    void substitute_in_edges(NodeId orig_to, NodeId new_to);

    /**
     * Get a non-owning view of the node with the given id. The view is
     * invalidated by any modification to the graph, including adding or removing nodes or
     * edges.
     *
     * @param id The id of the node to retrieve.
     * @return A view of the node with the given id.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    [[nodiscard]] NodeView node(NodeId id) const;

    /**
     * Get a non-owning reverse view of the node with the given id. The view is
     * invalidated by any modification to the graph, including adding or
     * removing nodes or edges. The reverse view has the same sequence as the
     * original node but in reverse orientation, and the incoming and outgoing
     * edges are swapped.
     *
     *
     * @param id The id of the node to retrieve.
     * @return A view of the node with the given id.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    [[nodiscard]] NodeView node_rev(NodeId id) const;

    /**
     * Get NodeIdRange corresponding to all currently nodes in the graph.
     *
     * @return A NodeIdRange of all nodes in the graph.
     */
    [[nodiscard]] NodeIdRange nodes() const;

    /**
     * Get the number of nodes in the graph.
     *
     * @return The total number of nodes in the graph.
     */
    [[nodiscard]] size_t nnodes() const;

    /**
     * Check if the node with the given id is a source node (i.e., has no incoming edges).
     *
     * @param id The id of the node to check.
     * @return true if the node is a source node, false otherwise.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    [[nodiscard]] bool is_source(NodeId id) const;

    /**
     * Check if the node with the given id is a sink node (i.e., has no outgoing edges).
     *
     * @param id The id of the node to check.
     * @return true if the node is a sink node, false otherwise.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    [[nodiscard]] bool is_sink(NodeId id) const;

    /**
     * @brief Return the length of the sequence associated with the given node.
     *
     * @param id The id of the node to check.
     * @return The length of the sequence associated with the given node.
     * @throws InvalidNodeException if the given id does not correspond to an
     * existing node.
     */
    [[nodiscard]] int node_size(NodeId id) const;

    /**
     * Get NodeIdRange corresponding to the source nodes in the graph. A source
     * node is a node with no incoming edges.
     *
     * @return A NodeIdRange of the source nodes in the graph.
     */
    [[nodiscard]] NodeIdRange source_nodes() const;

    /**
     * Get NodeIdRange corresponding to the sink nodes in the graph. A sink
     * node is a node with no outgoing edges.
     *
     * @return A NodeIdRange of the sink nodes in the graph.
     */
    [[nodiscard]] NodeIdRange sink_nodes() const;

    /**
     * Print the graph for debugging purposes.
     *
     */
    void print_graph();

private:
    std::unique_ptr<GraphImpl> impl_;
};

}  // namespace theseus
