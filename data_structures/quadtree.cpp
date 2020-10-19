// quadtree.cpp

#include <queue>

#include "quadtree.hpp"

int const MAX_LEAF_SIZE = 8;

using coord_t = spatial::coord_t;
using code_t = spatial::code_t;
using index_t = spatial::index_t;

template<typename T>
spatial::Quadtree<T>::Quadtree(coord_t x0, coord_t x1, coord_t y0, coord_t y1):
    root(new Node(NULL, 0, 0, (Rectangle){x0, x1, y0, y1})) {}

template<typename T>
spatial::Quadtree<T>::~Quadtree() { recursive_deconstruct(root); }

template<typename T>
spatial::Quadtree<T>::Node::Node(Node* p, int d, code_t c, Rectangle b): 
    depth(d), parent(p), code(c), bounds(b), center(midpoint(b)) {}

template<typename T>
spatial::Range spatial::Quadtree<T>::build(std::vector<T> const raw_data) {
    std::vector<Datum<T>> formatted_data;
    for (auto const& raw_datum : raw_data) {
        Point const new_point = {raw_datum[0], raw_datum[1]};
        Datum<T> const new_datum = {raw_datum, new_point};
        formatted_data.push_back(new_datum);
    }
    return recursive_build(root, formatted_data); 
}

template<typename T>
spatial::Range spatial::Quadtree<T>::recursive_build(
    Node* const node, std::vector<Datum<T>> data
) {
    if (data.size() <= MAX_LEAF_SIZE) {
        // Create a new leaf node
        index_t const index = leaves.size();
        node->leaf_range = {index, index};
        leaves.push_back((Leaf){node, data}); 
        return {index, index};
    } else {
        // Partition the data into four quadrants
        std::array<std::vector<Datum<T>>, 4> partition;
        for (auto const& datum : data) {
            int const quadrant = get_quadrant(node->center, datum.point);
            partition[quadrant].push_back(datum);
        }
        // Create 4 child-quadrants and recurse with the appropriate partition
        node->create_children();

        // The NW child will contain the leaf with the lowest Z-order code
        Range child_leaf_range;
        child_leaf_range = recursive_build(node->children[0], partition[0]);
        node->leaf_range.start = child_leaf_range.start;

        // The leaves in the NE and SW children all fall inside this range
        recursive_build(node->children[1], partition[1]);
        recursive_build(node->children[2], partition[2]);

        // The SE child will contain the leaf with the highest Z-order code
        child_leaf_range = recursive_build(node->children[3], partition[3]);
        node->leaf_range.end = child_leaf_range.end;

        return node->leaf_range;
    }
}

template<typename T>
void spatial::Quadtree<T>::recursive_deconstruct(Node* const node) {
    if (!node->is_leaf()) {
        for (auto const& child : node->children) {
            recursive_deconstruct(child);
        }
    }
    delete node;
}

/**
 * Query the quadtree for the k-nearest neighbours of a query point.
 * We do this with a doubly-greedy approach. First, we keep track of a priority
 * queue of quadtree nodes (initialized with just the root), and repeatedly
 * expand the node closest to the query point. As we find leaf nodes, we fill
 * up a second priority queue containing the nearest neighbours so far. Once
 * we have k points, and our "worst" nearest neighbour is closer than the next 
 * closest node, we terminate the search and return the k-nearest neighbours.
 */
template<typename T>
std::vector<T> spatial::Quadtree<T>::query_knn(
    unsigned const k, coord_t const x, coord_t const y
) {
    Point const query_point = {x, y};

    // Set up the node priority queue
    auto const node_cmp = [](NodePQE a, NodePQE b) {
        return (a.dist > b.dist);
    };
    std::priority_queue<NodePQE, std::vector<NodePQE>, decltype(node_cmp)> 
        node_pq(node_cmp);
    node_pq.push((NodePQE){root, distance(query_point, root->bounds)});

    // Set up the datum priority queue
    auto const datum_cmp = [](DatumPQE a, DatumPQE b) {
        return (a.dist < b.dist);
    };
    std::priority_queue<DatumPQE, std::vector<DatumPQE>, decltype(datum_cmp)> 
        datum_pq(datum_cmp);

    // Search the quadtree
    while (datum_pq.size() < k
           || datum_pq.top().dist > node_pq.top().dist) {
        NodePQE const npqe = node_pq.top();
        node_pq.pop();

        if (npqe.node->is_leaf()) {
            // Examine the points in the leaf
            index_t const leaf_index = npqe.node->leaf_range.start;
            for (auto const& datum : leaves[leaf_index].bucket) {
                coord_t const new_dist = distance(query_point, datum.point);
                // We have fewer than k points: take the new point
                if (datum_pq.size() < k) {
                    datum_pq.push((DatumPQE){datum, new_dist});
                }
                // The new point is better than our worst point: replace it
                else if (datum_pq.top().dist > new_dist) {
                    datum_pq.pop();
                    datum_pq.push((DatumPQE){datum, new_dist});
                }
            }
        } else {
            // Expand the node, adding its children to the priority queue
            for (auto const& child : npqe.node->children) {
                coord_t const dist = distance(query_point, child->bounds);
                node_pq.push((NodePQE){child, dist});
            }
        }
    }
    std::vector<T> query_bucket;
    query_bucket.reserve(k);
    while (!datum_pq.empty()) {
        DatumPQE e = datum_pq.top();
        query_bucket.push_back(e.datum.data);
        datum_pq.pop();
    }
    return query_bucket;
}

/**
 * A note on the differing comparison operators: Z-ordering 0 is to the NW,
 * but the geographic coordinate (0,0) is in the SW. So, while it looks weird,
 * this is just converting from geographic coordinates to Z-ordering.
 */
template<typename T>
int spatial::Quadtree<T>::get_quadrant(Point origin, Point p) {
    return ((p.x > origin.x) + ((p.y < origin.y) << 1));
}

template<typename T>
int spatial::Quadtree<T>::num_leaves() { return leaves.size(); }

/**
 * Verify that EVERY leaf in the quadtree is at depth k
 * Note that this will only pass for a complete quadtree of depth k
 */
template<typename T>
bool spatial::Quadtree<T>::check_depth_equals(unsigned const k) {
    for (auto const& leaf : leaves) {
        if (leaf.node->depth != k) return false;
    }
    return true;
}

template<typename T>
void spatial::Quadtree<T>::Node::create_children() {
    Rectangle NW_bounds = {bounds.xmin, center.x, center.y, bounds.ymax};
    children[0] = new Node(this, depth+1, (code << 2) + 0, NW_bounds);

    Rectangle NE_bounds = {center.x, bounds.xmax, center.y, bounds.ymax};
    children[1] = new Node(this, depth+1, (code << 2) + 1, NE_bounds);

    Rectangle SW_bounds = {bounds.xmin, center.x, bounds.ymin, center.y};
    children[2] = new Node(this, depth+1, (code << 2) + 2, SW_bounds);
    
    Rectangle SE_bounds = {center.x, bounds.xmax, bounds.ymin, center.y};
    children[3] = new Node(this, depth+1, (code << 2) + 3, SE_bounds);
}

template<typename T>
bool spatial::Quadtree<T>::Node::is_leaf() { 
    return (leaf_range.start == leaf_range.end);
}
