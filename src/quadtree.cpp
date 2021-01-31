// quadtree.cpp

#include <queue>

#include "quadtree.hpp"

/**
 * These two constants are mutually exclusive:
 *  LEAF_CAPACITY is used for insert-based construction
 *  TARGET_DEPTH is used for bulk loading
 */
int const LEAF_CAPACITY = 16;

using coord_t = spatial::coord_t;
using code_t = spatial::code_t;
using index_t = spatial::index_t;

/**
 * We nudge the maximum bounds of the quadtree by a tiny amount to fix a
 * rare edge case in zorder_hash() where a point is right on the boundary.
 * (a point on the boundary results in a Z-order code outside of the tree)
 */
template<typename T>
spatial::Quadtree<T>::Quadtree(coord_t x0, coord_t x1, coord_t y0, coord_t y1):
    root(std::make_unique<Node>(0, 0, (Rectangle){x0, x1+0.01, y0, y1+0.01}))
{ }

template<typename T>
spatial::Quadtree<T>::Node::Node(int d, code_t c, Rectangle b): 
    depth(d), 
    code(c), 
    bounds(b), 
    center(midpoint(b)) 
{ }

/**
 * Now that we have two methods with which to build the quadtree,
 * we'll use build() as an alias for bulk load, for compatibiliy's sake
 */
template<typename T>
void spatial::Quadtree<T>::build(std::vector<T> const& raw_data) {
    root->insert(*this, datumize<T>(raw_data));
}

/**
 * Recursively insert a collection of data into the quadtree.
 */
template<typename T>
spatial::Range spatial::Quadtree<T>::Node::insert(
    Quadtree<T>& tree, std::vector<Datum<T>> data
) {
    if (data.size() <= LEAF_CAPACITY) {
        // Create a new leaf node
        index_t const leaf_idx = tree.leaves.size();
        this->leaf_range = {leaf_idx, leaf_idx};
        tree.leaves.push_back(data); 
        return {leaf_idx, leaf_idx};
    } else {
        // Partition the data into four quadrants
        std::array<std::vector<Datum<T>>, 4> partition;
        for (auto const& datum : data) {
            int const quadrant = get_quadrant(datum.point);
            partition[quadrant].push_back(datum);
        }
        // Create 4 child-quadrants and recurse with the appropriate partition
        this->create_children();

        // The NW child will contain the leaf with the lowest Z-order code
        Range child_leaf_range = children[0]->insert(tree, partition[0]);
        this->leaf_range.start = child_leaf_range.start;

        // The leaves in the NE and SW children all fall inside this range
        children[1]->insert(tree, partition[1]);
        children[2]->insert(tree, partition[2]);

        // The SE child will contain the leaf with the highest Z-order code
        child_leaf_range = children[3]->insert(tree, partition[3]);
        this->leaf_range.end = child_leaf_range.end;

        return this->leaf_range;
    }
}

/**
 * k-nearest neighbour query using our distance browsing algorithm
 */
template<typename T>
std::vector<T> spatial::Quadtree<T>::query_knn(
    unsigned const k, coord_t const x, coord_t const y
) const {
    Point const query_point = {x, y};

    NodePQ node_pq(query_point);
    node_pq.push(root.get());
    DatumPQ datum_pq(query_point);

    while (!node_pq.empty() && (
        datum_pq.size() < k
        || datum_pq.peek().dist > node_pq.peek().dist
    )) {
        Node* next_node = node_pq.pop().node;
        if (next_node->is_leaf()) {
            for (auto const& datum : leaves[next_node->leaf_range.start]) {
                if (datum_pq.size() < k) {
                    datum_pq.push(datum);
                }
                else {
                    datum_pq.choose(datum);
                }
            }
        } else {
            node_pq.expand(next_node);
        }
    } 

    std::vector<T> query_bucket;
    query_bucket.reserve(k);
    while (!datum_pq.empty()) {
        query_bucket.push_back(datum_pq.pop().datum.data);
    }
    return query_bucket;
}

/**
 * A note on the differing comparison operators: Z-ordering 0 is to the NW,
 * but the geographic coordinate (0,0) is in the SW. So, while it looks weird,
 * this is just converting from geographic coordinates to Z-ordering.
 */
template<typename T>
int spatial::Quadtree<T>::Node::get_quadrant(Point const p) const {
    return ((p.x > center.x) + ((p.y > center.y) << 1));
}

template<typename T>
int spatial::Quadtree<T>::num_leaves() const { return leaves.size(); }

template<typename T>
void spatial::Quadtree<T>::Node::create_children() {
    Rectangle const SW_bounds = {bounds.xmin, center.x, bounds.ymin, center.y};
    children[0] = std::make_unique<Node>(depth+1, (code << 2) + 0, SW_bounds);

    Rectangle const SE_bounds = {center.x, bounds.xmax, bounds.ymin, center.y};
    children[1] = std::make_unique<Node>(depth+1, (code << 2) + 1, SE_bounds);

    Rectangle const NW_bounds = {bounds.xmin, center.x, center.y, bounds.ymax};
    children[2] = std::make_unique<Node>(depth+1, (code << 2) + 2, NW_bounds);

    Rectangle const NE_bounds = {center.x, bounds.xmax, center.y, bounds.ymax};
    children[3] = std::make_unique<Node>(depth+1, (code << 2) + 3, NE_bounds);
}

template<typename T>
bool spatial::Quadtree<T>::Node::is_leaf() const { 
    return (leaf_range.start == leaf_range.end);
}
