// quadtree.cpp

#include "quadtree.hpp"

using namespace std;

const int MAX_LEAF_SIZE = 8;

template<typename T>
spatial::Quadtree<T>::Quadtree(double x0, double x1, double y0, double y1):
    root(new Node(NULL, 0, 0, (Rectangle){x0, x1, y0, y1})) {}

template<typename T>
spatial::Quadtree<T>::Node::Node(Node* p, int d, long long int c, Rectangle b): 
    parent(p), depth(d), code(c), bounds(b), center(midpoint(b)) {}

template<typename T>
spatial::Range spatial::Quadtree<T>::build(vector<T> raw_data) {
    vector<Datum<T>> formatted_data;
    for (int i=0; i<(int)raw_data.size(); i++) {
        Point new_point = {raw_data[i][0], raw_data[i][1]};
        Datum<T> new_datum = {raw_data[i], new_point};
        formatted_data.push_back(new_datum);
    }
    return recursive_build(root, formatted_data); 
}

template<typename T>
spatial::Range spatial::Quadtree<T>::recursive_build(
        Node* node, vector<Datum<T>> data) {
    if (data.size() <= MAX_LEAF_SIZE) {
        // Create a new leaf node
        int index = leaves.size();
        node->leaf_range = {index, index};
        Leaf new_leaf = {node, data};
        leaves.push_back(new_leaf); 
        return {index, index};
    } else {
        // Partition the data
        array<vector<Datum<T>>, 4> partition;
        for (int i=0; i<(int)data.size(); i++) {
            int quadrant = get_quadrant(node->center, data[i].point);
            partition[quadrant].push_back(data[i]);
        }
        data.clear(); // Free up some memory

        // Create 4 child-quadrants and recurse with the appropriate partition
        node->create_children();
        Range child_leaf_range;

        // The NW child will contain the leaf with the lowest Z-order code
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

/**
 * A note on the differing comparison operators: Z-ordering 0 is to the NW,
 * but the geographic coordinate (0,0) is in the SW. So, while it looks weird,
 * this is just us converting from geographic coords to Z-ordering.
 */
template<typename T>
int spatial::Quadtree<T>::get_quadrant(Point origin, Point p) {
    return ((p.x > origin.x) + ((p.y < origin.y) << 1));
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
int spatial::Quadtree<T>::num_leaves() { return leaves.size(); }
