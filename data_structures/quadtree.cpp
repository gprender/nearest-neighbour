// quadtree.cpp

#include "quadtree.hpp"

using namespace std;

const int MAX_LEAF_SIZE = 8;

template<typename T>
spatial::Quadtree<T>::Quadtree(double x0, double x1, double y0, double y1) {
    height = 1;
    load = 0;
    bounds = (spatial::Rectangle){x0, x1, y0, y1};
    xmedian = (bounds.xmin + bounds.xmax) / 2;
    ymedian = (bounds.ymin + bounds.ymax) / 2;
    try { 
        data = new spatial::Datum<T>[MAX_LEAF_SIZE];
    } catch (...) {
        cout << "Exception thrown while allocating a new datum\n";
        exit(0);
    }
    NW = NE = SW = SE = NULL;
}

template<typename T>
int spatial::Quadtree<T>::insert(T data) {
    struct spatial::Point new_point = (spatial::Point){data[0], data[1]};
    struct spatial::Datum<T> new_datum = (spatial::Datum<T>){data, new_point};
    return insert(new_datum);
}

template<typename T>
int spatial::Quadtree<T>::insert(spatial::Datum<T> datum) {
    if (isLeaf()) {
        this->data[load] = datum;
        load++;
        if (getLoad() >= MAX_LEAF_SIZE) height = split_node();
    } else { 
        int subtree_height = recursive_insert(datum);
        if (subtree_height >= height) height = 1 + subtree_height;
        load++;
    }
    return height;
}

template<typename T>
int spatial::Quadtree<T>::recursive_insert(spatial::Datum<T> datum) {
    double x = datum.point.x;
    double y = datum.point.y;
    int subtree_height;

    // Determine which quadrant to recurse on
    if (x < xmedian) { 
        if (y < ymedian) {
            subtree_height = SW->insert(datum);
        } else {
            subtree_height = NW->insert(datum);
        }
    } else {
        if (y < ymedian) {
            subtree_height = SE->insert(datum);
        } else {
            subtree_height = NE->insert(datum);
        }
    }
    /**
     * We could add a 3rd case here to check that the datum
     * actually falls inside of the current cell's bounds
     */
    return subtree_height;
}

// Split the current full node into 4 quadrants
template<typename T>
int spatial::Quadtree<T>::split_node() {
    NW = new spatial::Quadtree<T>(bounds.xmin, xmedian, ymedian, bounds.ymax);
    NE = new spatial::Quadtree<T>(xmedian, bounds.xmax, ymedian, bounds.ymax);
    SW = new spatial::Quadtree<T>(bounds.xmin, xmedian, bounds.ymin, ymedian);
    SE = new spatial::Quadtree<T>(xmedian, bounds.xmax, bounds.ymin, ymedian);

    // Distribute the node's points to the appropriate children
    int max_subtree_height = 1;
    int temp_subtree_height;
    spatial::Datum<T> datum;
    for (int i=0; i<MAX_LEAF_SIZE; i++) {
        datum = this->data[i];
        temp_subtree_height = recursive_insert(datum);
        if (temp_subtree_height > max_subtree_height) {
            max_subtree_height = temp_subtree_height;
        }
    }
    delete[] this->data;
    return (1 + max_subtree_height);
}

template<typename T>
int spatial::Quadtree<T>::getHeight() { return height; }

template<typename T>
int spatial::Quadtree<T>::getLoad() { return load; }

template<typename T>
bool spatial::Quadtree<T>::isLeaf() { return (NW == NULL); }
