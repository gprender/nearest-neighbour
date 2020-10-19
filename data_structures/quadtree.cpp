// quadtree.cpp

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
        Node* const node, std::vector<Datum<T>> data) {
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
 * Query the quadtree for all points "around" a given query point.
 * Currently, this means all points in the query point's leaf, as well as
 * all points in the 8 leaves which border on the query point's leaf.
 */
template<typename T>
std::vector<T> spatial::Quadtree<T>::query(coord_t const x, coord_t const y) {
    Node* const origin_node = find_leaf((Point){x, y});
    std::vector<code_t> neighbour_codes = calc_neighbour_codes(
        origin_node->code, origin_node->depth);

    // Start by collecting points from the origin node and its siblings
    // No traversal needed, since all four are contiguous in the leaf vector
    std::vector<T> query_bucket;
    index_t const origin_index = origin_node->leaf_range.start;
    int const origin_quadrant = origin_node->code & 3;
    for (int i=(origin_index - origin_quadrant); i<4; i++) {
        for (auto const& datum : leaves[i].bucket) {
            query_bucket.push_back(datum.data);
        }
    }
    // Next, collect points from the five non-sibling neighbours
    // A small optimization here could be to group sibling codes together
    Node* last_node = NULL;
    for (auto const target_code : neighbour_codes) {
        Node* const target_node = traverse(origin_node, target_code);

        // Ensure we're not querying a node >1 time; this can happen if we
        // reach a leaf before getting to the target depth
        if (target_node != last_node) {
            // Since we keep track of each nodes' contiguous range of leaves
            // in the leaf vector, collecting points from a leaf is the same
            // as collecting points from an internal node
            index_t const start = target_node->leaf_range.start;
            index_t const end = target_node->leaf_range.end;
            for (index_t i=start; i<=end; i++) {
                for (auto const& datum : leaves[i].bucket) {
                    query_bucket.push_back(datum.data);
                }
            }
            last_node = target_node;
        }
    }
    return query_bucket;
}

/**
 * Calculate the neighbouring location codes of a given origin location code.
 * This function uses a finite state machine (defined in the header file) to
 * calculate neighbour codes for each direction. The algorithm is described
 * in detail in by Yoder and Bloniarz (2006).
 */
template<typename T>
std::vector<code_t> spatial::Quadtree<T>::calc_neighbour_codes(
        const code_t origin_code, const int depth) {
    // Calculate a code for each direction, clockwise starting from north
    std::vector<code_t> codes;
    for (int i=0; i<8; i++) {
        int direction = i; // Initial neighbour direction
        code_t c0 = origin_code;

        // Calculate each new digit in order (right to left)
        for (int j=0; j<depth; j++) {
            code_t const d0 = (origin_code >> (2*j)) & 3; // origin digit j
            code_t const d1 = fsm[d0][direction][0]; // neighbour digit j
            code_t const c1 = ~(3 << (2*j)); // bitmask to clear digit j

            c0 &= c1; // Clear the jth digit of the origin code
            c0 |= (d1 << (2*j)); // Set the jth digit of the new neighbour code

            if (fsm[d0][direction][1] == -1){
                break; // We've found a sibling -> halt
            } else {
                direction = fsm[d0][direction][1]; // Update travel direction
            }
            // If we reach the top level of the quadtree without finding
            // a sibling, then a neighbour doesn't exist in this direction
            if (j == depth-1) c0 = -1;
        }
        // Conditions for filtering out certain types of neighbour codes
        if (c0 == -1) continue; // Neighbour does not exist
        if (c0 >> 2 == origin_code >> 2) continue; // Neighbour is a sibling

        codes.push_back(c0);
    }
    return codes;
}

// Traverse the quadtree from an origin node to a target location code.
template<typename T>
typename spatial::Quadtree<T>::Node* spatial::Quadtree<T>::traverse(
        Node* const origin_node, code_t const target_code) {
    int const origin_depth = origin_node->depth;
    
    // Ascend the tree until we find a common ancestor node
    Node* current_node = origin_node;
    while (true) {
        int const depth_difference = origin_depth - current_node->depth;
        if (current_node->code == target_code >> (2*depth_difference)) {
            break;
        } else {
            current_node = current_node->parent;
        }
    }

    // Descend the tree according to the target location code
    // Note that we might hit a leaf before reaching the target location,
    // or we might discover that the target is an internal node
    while (!current_node->is_leaf()) {
        int const depth_difference = origin_depth - current_node->depth;
        int const next_quadrant = (target_code >> (2*(depth_difference - 1))) & 3;
        current_node = current_node->children[next_quadrant];
    }
    return current_node;
}

// Find the leaf whose bounds contain the given point
template<typename T>
typename spatial::Quadtree<T>::Node* spatial::Quadtree<T>::find_leaf(const Point p) {
    Node* current_node = root;
    int next_quadrant;
    while (!current_node->is_leaf()) {
        next_quadrant = get_quadrant(current_node->center, p);
        current_node = current_node->children[next_quadrant];
    }
    return current_node;
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
