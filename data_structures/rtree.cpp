// rtree.cpp

#include <iostream>
#include <memory>

#include "rtree.hpp"

using coord_t = spatial::coord_t;
using code_t = spatial::code_t;
using index_t = spatial::index_t;

int const M = 8;

template<typename T>
spatial::Rtree<T>::Rtree(): 
    root_entry(std::make_shared<InternalEntry>()) 
{
    root_entry->bounding_box = (Rectangle){0, 0, 0, 0};
    root_entry->node = std::make_shared<Node>(*this);
}

// ~Rtree() {}

template<typename T>
spatial::Rtree<T>::Node::Node(Rtree& rt): 
    load(0), 
    m(0), 
    rtree(rt) 
{ }

// ~Node() {}

/**
 * Construct an R-tree from the given point data.
 * Currently, this is just doing point-by-point insertion
 */
template<typename T>
void spatial::Rtree<T>::build(std::vector<T> const raw_data) {
    // Pick an inital bounding box for the root
    Rectangle initial_seed = {
        raw_data[0][0], raw_data[0][0],
        raw_data[0][1], raw_data[0][1]
    };
    root_entry->bounding_box = initial_seed;

    // Insert each data point into the R-tree
    data.reserve(raw_data.size());
    for (auto const& raw_datum : raw_data) {
        Point const new_point = {raw_datum[0], raw_datum[1]};
        Datum<T> const new_datum = {raw_datum, new_point};
        insert(new_datum);
    } 
}

template<typename T>
void spatial::Rtree<T>::insert(Datum<T> new_datum) {
    data.push_back(new_datum);

    // Expand the root's bounding box if necessary
    root_entry->bounding_box = min_bounding_box(
        root_entry->bounding_box, new_datum.point
    );

    // Recursively insert the point into the root node
    if (root_entry->node->insert(data.size()-1)) {
        split_root();
    }
}

/**
 * Recursively insert a point into the current node.
 * If a node exceeds M entries, we return 'true' to indicate that a split
 * is required, since splitting happens at the parent's level.  
 */
template<typename T>
bool spatial::Rtree<T>::Node::insert(index_t const point_idx) {
    Point const p = rtree.data[point_idx].point;
    if (this->is_leaf()) {
        // add the point to the current node
        auto new_entry = std::make_shared<LeafEntry>();
        new_entry->bounding_box = (Rectangle){ p.x, p.x, p.y, p.y };
        new_entry->idx = point_idx;
        entries.push_back(new_entry);
        m++;
    } else {
        // we're in an internal node, and need to descend further
        int const branch_idx = choose_branch(p);
        auto entry = std::dynamic_pointer_cast<InternalEntry>(
            entries[branch_idx]
        );

        // expand the child's bounding box as needed
        entry->bounding_box = min_bounding_box(
            entry->bounding_box, p
        );

        // recurse on the child node
        if (entry->node->insert(point_idx)) {
            split(branch_idx);  // split if the child overflows
        } 
    }
    load++;
    return (m > M);  // if (m > M), this branch needs to be split
}

/**
 * Split an overflowing entry (aka child node) into two new entries.
 * Note that the object we're calling split() on is the PARENT of the node
 * to be split, not the node itself.
 */
template<typename T>
void spatial::Rtree<T>::Node::split(int const branch_idx) {
    // pop the overflowing branch from the 'entries' vector

    // pick some seed entries using the split heuristic

    // push the newly seeded branches onto 'entries'

    // distribute the leftover points
}

/**
 * When the root node overflows, we need some special logic, 
 * since it has no parent node.
 */
template<typename T>
void spatial::Rtree<T>::split_root() {
    // make a new root, with the old root being its only entry

    // root->split(0);
}

/**
 * Pick an appropriate child node for a point.
 * Specifically, we pick the bounding box which requires the
 * smallest (area) expansion to accommodate the new point.
 */
template<typename T>
int spatial::Rtree<T>::Node::choose_branch(Point const p) const {
    coord_t min_expansion;
    int best_choice = -1;
    int pos = 0;
    for (auto& entry : entries) {
        Rectangle const expanded_bb = min_bounding_box(
            entry->bounding_box, p
        );
        coord_t const current_expansion = (
            area(expanded_bb) - area(entry->bounding_box)
        );
        if (best_choice == -1 || current_expansion < min_expansion) {
            min_expansion = current_expansion;
            best_choice = pos;
        }
        if (min_expansion == 0) break; // We've found the best case
        pos++;
    } 
    return best_choice;
}

template<typename T>
spatial::index_t spatial::Rtree<T>::get_load() const { 
    return root_entry->node->load; 
}

template<typename T>
bool spatial::Rtree<T>::Node::is_leaf() const { 
    return load < M; 
}
