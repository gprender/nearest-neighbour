// quadtree.cpp

#include <queue>

#include "quadtree.hpp"

int const LEAF_CAPACITY = 64;

using coord_t = spatial::coord_t;
using code_t = spatial::code_t;
using index_t = spatial::index_t;

template<typename T>
spatial::Quadtree<T>::Quadtree(coord_t x0, coord_t x1, coord_t y0, coord_t y1):
    root(std::make_unique<Node>(0, 0, (Rectangle){x0, x1, y0, y1}))
{ }

template<typename T>
spatial::Quadtree<T>::~Quadtree() 
{ }

template<typename T>
spatial::Quadtree<T>::Node::Node(int d, code_t c, Rectangle b): 
    depth(d), 
    code(c), 
    bounds(b), 
    center(midpoint(b)) 
{ }

template<typename T>
spatial::Quadtree<T>::Node::~Node() 
{ }

template<typename T>
spatial::Range spatial::Quadtree<T>::build(std::vector<T> const& raw_data) {
    std::vector<Datum<T>> formatted_data;
    for (auto const& raw_datum : raw_data) {
        Point const new_point = {raw_datum[0], raw_datum[1]};
        Datum<T> const new_datum = {raw_datum, new_point};
        formatted_data.push_back(new_datum);
    }
    return root->insert(*this, formatted_data);
}

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

    while (
        datum_pq.size() < k
        || datum_pq.peek().dist > node_pq.peek().dist
    ) {
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
    return ((p.x > center.x) + ((p.y < center.y) << 1));
}

template<typename T>
int spatial::Quadtree<T>::num_leaves() const { return leaves.size(); }

/**
 * Verify that EVERY leaf in the quadtree is at depth k
 * Note that this will only pass for a complete quadtree of depth k
 */
template<typename T>
bool spatial::Quadtree<T>::depth_equals(unsigned const k) const {
    for (auto const& leaf : leaves) {
        if (leaf.node->depth != k) return false;
    }
    return true;
}

template<typename T>
void spatial::Quadtree<T>::Node::create_children() {
    Rectangle const NW_bounds = {bounds.xmin, center.x, center.y, bounds.ymax};
    children[0] = std::make_unique<Node>(depth+1, (code << 2) + 0, NW_bounds);

    Rectangle const NE_bounds = {center.x, bounds.xmax, center.y, bounds.ymax};
    children[1] = std::make_unique<Node>(depth+1, (code << 2) + 1, NE_bounds);

    Rectangle const SW_bounds = {bounds.xmin, center.x, bounds.ymin, center.y};
    children[2] = std::make_unique<Node>(depth+1, (code << 2) + 2, SW_bounds);
    
    Rectangle const SE_bounds = {center.x, bounds.xmax, bounds.ymin, center.y};
    children[3] = std::make_unique<Node>(depth+1, (code << 2) + 3, SE_bounds);
}

template<typename T>
bool spatial::Quadtree<T>::Node::is_leaf() const { 
    return (leaf_range.start == leaf_range.end);
}
