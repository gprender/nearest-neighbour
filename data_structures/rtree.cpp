// rtree.cpp

#include <iostream>
#include <memory>

#include "rtree.hpp"

using coord_t = spatial::coord_t;
using area_t = spatial::area_t;
using index_t = spatial::index_t;

int const M = 8;

template<typename T>
spatial::Rtree<T>::Rtree():
    root_entry(std::make_unique<Entry>(
        (Rectangle){0,0,0,0}, 
        std::make_shared<Node>()
    ))
{ }

template<typename T>
spatial::Rtree<T>::~Rtree() { }

template<typename T>
spatial::Rtree<T>::Node::Node(): 
    load(0)
{ }

template<typename T>
spatial::Rtree<T>::Node::~Node() { }

/**
 * Construct an R-tree from the given point data.
 * Currently, this is just doing point-by-point insertion
 */
template<typename T>
void spatial::Rtree<T>::build(std::vector<T> const& raw_data) {
    // Pick an inital bounding box for the root
    Rectangle initial_seed = {
        raw_data[0][0], raw_data[0][0],
        raw_data[0][1], raw_data[0][1]
    };
    root_entry->set_mbb(initial_seed);

    // Insert each data point into the R-tree
    data.reserve(raw_data.size());
    for (auto const& raw_datum : raw_data) {
        Point const new_point = {raw_datum[0], raw_datum[1]};
        Datum<T> const new_datum = {raw_datum, new_point};
        insert(new_datum);
    } 
}

template<typename T>
void spatial::Rtree<T>::insert(Datum<T> const& new_datum) {
    data.push_back(new_datum);

    // Expand the root's bounding box if necessary
    root_entry->set_mbb(min_bounding_box(
        root_entry->get_mbb(), new_datum.point
    ));

    // Recursively insert the point into the root node
    if (root_entry->get_node()->insert(data.back())) {
        split_root();  // split the root if it overflows
    }
}

/**
 * Recursively insert a point into the current node.
 * If a node exceeds M entries, we return 'true' to indicate that a split
 * is required, since splitting happens at the parent's level.  
 */
template<typename T>
bool spatial::Rtree<T>::Node::insert(Datum<T> const& datum) {
    Point const p = datum.point;
    if (this->is_leaf()) {
        // add the point to the current node
        entries.push_back(Entry(
            (Rectangle){p.x, p.x, p.y, p.y},
            std::make_shared<Datum<T>>(datum)
        ));
    } else {
        // we're in an internal node, and need to descend further
        int const branch_idx = choose_branch(p);
        Entry child_entry = entries[branch_idx];

        // expand the child's bounding box as needed
        child_entry.set_mbb(min_bounding_box(
            child_entry.get_mbb(), p
        ));

        // recurse on the child node
        if (child_entry.get_node()->insert(datum)) {
            split(branch_idx);  // split if the child overflows
        } 
    }
    load++;
    return (entries.size() > M);  // check for node overflow
}

/**
 * When the root node overflows, we need some special logic, 
 * since it has no parent node.
 */
template<typename T>
void spatial::Rtree<T>::split_root() {
    // make a new root, with the old root being its only entry
    auto other_entry = std::make_unique<Entry>(
        root_entry->get_mbb(),
        std::make_shared<Node>()
    );
    root_entry.swap(other_entry);
    // other_entry now contains the OLD root entry
    root_entry->get_node()->entries.push_back(*other_entry);
    root_entry->get_node()->load = other_entry->get_node()->load;
    // root node entries[0] is now the old, overflowing root
    root_entry->get_node()->split(0);
}

/**
 * Split an overflowing entry (aka child node) into two new entries.
 * Note that the object we're calling split() on is the PARENT of the node
 * to be split, not the node itself.
 */
template<typename T>
void spatial::Rtree<T>::Node::split(int const branch_idx) {
    // pop the overflowing branch from the 'entries' vector
    auto const overflowing_node = entries[branch_idx].get_node();
    entries.erase(entries.begin() + branch_idx);

    // make some seed entries using the split heuristic
    pick_seeds(overflowing_node->entries);

    // distribute the leftover child entries between the seeds
    distribute(overflowing_node->entries);
}

/**
 * Pick a number of "good" seed MBBs from the given choices.
 * Currently, this function uses the quadratic split heuristic.
 */
template<typename T>
void spatial::Rtree<T>::Node::pick_seeds(
    std::vector<Entry> const& entry_choices
) {
    unsigned best_e1 = -1, best_e2 = -1;
    area_t max_d = -1;
    for (unsigned i=0; i<entry_choices.size(); i++) {
        Entry const e1 = entry_choices[i];
        for (unsigned j=i+1; j<entry_choices.size(); j++) {
            Entry const e2 = entry_choices[j];
            Rectangle const r = min_bounding_box(e1.get_mbb(), e2.get_mbb());
            area_t const d = area(r) - area(e1.get_mbb()) - area(e2.get_mbb());
            if (d > max_d) {
                max_d = d;
                best_e1 = i;
                best_e2 = j;
            }
        }
    }
    // create two new entries corresponding to the above seed MBBs
    entries.push_back(Entry(
        entry_choices[best_e1].get_mbb(),
        std::make_shared<Node>()
    ));
    entries.push_back(Entry(
        entry_choices[best_e2].get_mbb(),
        std::make_shared<Node>()
    ));
}

/**
 * Distribute leftover entries after splitting an overflowing node.
 * This function is a bit of a mouthful, particularly the if/else stack
 * at the end, so it could probably do with some revisiting.
 */
template<typename T>
void spatial::Rtree<T>::Node::distribute(
    std::vector<Entry>& leftover_entries
) {
    // our two "groups" are the child nodes that were just created
    Entry& g1 = entries[entries.size()-1];
    Entry& g2 = entries[entries.size()-2];
    while (!leftover_entries.empty()) {
        // pop the "best" entry from the leftovers vector
        int const next_idx = pick_next(leftover_entries);
        Entry const next_entry = leftover_entries[next_idx];
        leftover_entries.erase(leftover_entries.begin() + next_idx);

        // calculate the MBB expansion needed by each group
        Rectangle const g1_expanded_mbb = min_bounding_box(
            g1.get_mbb(), next_entry.get_mbb()
        );
        Rectangle const g2_expanded_mbb = min_bounding_box(
            g2.get_mbb(), next_entry.get_mbb()
        );

        area_t const g1_expansion = area(g1_expanded_mbb) - area(g1.get_mbb());
        area_t const g2_expansion = area(g2_expanded_mbb) - area(g2.get_mbb());

        // choose the group which requires the least expansion
        if (g1_expansion < g2_expansion) {
            g1.set_mbb(g1_expanded_mbb);
            g1.get_node()->entries.push_back(next_entry);
            g1.get_node()->load++;
        } else if (g1_expansion == g2_expansion) {
            // break the tie (pick the group w/ the smallest MBB)
            if (area(g1.get_mbb()) < area(g2.get_mbb())) {
                g1.set_mbb(g1_expanded_mbb);
                g1.get_node()->entries.push_back(next_entry);
                g1.get_node()->load++;
            } else {
                g2.set_mbb(g2_expanded_mbb);
                g2.get_node()->entries.push_back(next_entry);
                g2.get_node()->load++;
            }
        } else {
            g2.set_mbb(g2_expanded_mbb);
            g2.get_node()->entries.push_back(next_entry);
            g2.get_node()->load++;
        }
    }
}

/**
 * Pick the "best" leftover entry to distribute next.
 */
template<typename T>
int spatial::Rtree<T>::Node::pick_next(
    std::vector<Entry> const& leftover_entries
) const {
    Entry const& g1 = entries[entries.size()-1];
    Entry const& g2 = entries[entries.size()-2];

    area_t max_diff = 0;
    int pos = 0, best_choice = 0;
    for (auto const& entry : leftover_entries) {
        area_t const d1 = area(
            min_bounding_box(g1.get_mbb(), entry.get_mbb())      
        ) - area(g1.get_mbb());

        area_t const d2 = area(
            min_bounding_box(g2.get_mbb(), entry.get_mbb())      
        ) - area(g2.get_mbb());

        if (std::abs(d1 - d2) > max_diff) {
            max_diff = std::abs(d1 - d2);
            best_choice = pos;
        }
        pos++;
    }
    return best_choice;
}

/**
 * Pick an appropriate child node for a point.
 * Specifically, we pick the bounding box which requires the
 * smallest (area) expansion to accommodate the new point.
 */
template<typename T>
int spatial::Rtree<T>::Node::choose_branch(Point const p) const {
    area_t min_expansion = -1;
    int pos = 0, best_choice = -1;
    for (auto const& entry : entries) {
        Rectangle const expanded_bb = min_bounding_box(
            entry.get_mbb(), p
        );
        coord_t const current_expansion = (
            area(expanded_bb) - area(entry.get_mbb())
        );
        if (best_choice == -1 || current_expansion < min_expansion) {
            min_expansion = current_expansion;
            best_choice = pos;
        } 
        else if (current_expansion == min_expansion) {
            // break the tie (choose the MBB w/ the smallest area)
            if (area(entry.get_mbb()) < area(entries[best_choice].get_mbb())) {
                min_expansion = current_expansion;
                best_choice = pos;
            }
        }
        pos++;
    } 
    return best_choice;
}

template<typename T>
spatial::index_t spatial::Rtree<T>::get_load() const { 
    return root_entry->get_node()->load; 
}

template<typename T>
bool spatial::Rtree<T>::Node::is_leaf() const { 
    return (entries.size() == load);
}
