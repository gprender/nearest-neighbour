// zgrid.cpp

#include "zgrid.hpp"

using coord_t = spatial::coord_t;
using code_t = spatial::code_t;
using index_t = spatial::index_t;

template<typename T>
spatial::Zgrid<T>::Zgrid(coord_t x0, coord_t x1, coord_t y0, coord_t y1):
    root(std::make_unique<Node>(0, 0, (Rectangle){x0, x1+0.01, y0, y1+0.01}))
{ }

template<typename T>
spatial::Zgrid<T>::Node::Node(code_t c, int d, Rectangle b):
    code(c),
    depth(d), 
    bounds(b), 
    center(midpoint(b))
{ }

template<typename T>
void spatial::Zgrid<T>::build(std::vector<T> const& raw_data, int const r) {
    zgrid_bin(datumize<T>(raw_data), r);
}

template<typename T>
void spatial::Zgrid<T>::zgrid_bin(std::vector<Datum<T>> const& data, int const r) {
    grid.resize(std::pow(4,r));
    for (auto const& datum : data) {
        code_t zorder_code = zorder_hash(datum.point, r);
        grid[zorder_code].push_back(datum);
    }
    root->populate(r);
}

template<typename T>
std::vector<T> spatial::Zgrid<T>::query_knn(
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
            for (auto const& datum : grid[next_node->code]) {
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

template<typename T>
code_t spatial::Zgrid<T>::zorder_hash(Point const p, int const r) const {
    Rectangle const& b = root->bounds;
    int cellx = grid_index(p.x, b.xmin, b.xmax, std::pow(2,r));
    int celly = grid_index(p.y, b.ymin, b.ymax, std::pow(2,r));
    return interleave(cellx, celly);
}

template<typename T>
void spatial::Zgrid<T>::Node::populate(int const r) {
    if (r > 0) {
        create_children();
        for (int i=0; i<4; i++) { 
            children[i]->populate(r-1);
        }
    }
}

template<typename T>
void spatial::Zgrid<T>::Node::create_children() {
    Rectangle const SW_bounds = {bounds.xmin, center.x, bounds.ymin, center.y};
    children[0] = std::make_unique<Node>((code << 2) + 0, depth+1, SW_bounds);

    Rectangle const SE_bounds = {center.x, bounds.xmax, bounds.ymin, center.y};
    children[1] = std::make_unique<Node>((code << 2) + 1, depth+1, SE_bounds);

    Rectangle const NW_bounds = {bounds.xmin, center.x, center.y, bounds.ymax};
    children[2] = std::make_unique<Node>((code << 2) + 2, depth+1, NW_bounds);

    Rectangle const NE_bounds = {center.x, bounds.xmax, center.y, bounds.ymax};
    children[3] = std::make_unique<Node>((code << 2) + 3, depth+1, NE_bounds);
}

template<typename T>
bool spatial::Zgrid<T>::Node::is_leaf() const {
    if (children[0]) {
        return false;
    } else {
        return true;
    }
}

template<typename T>
size_t spatial::Zgrid<T>::size() {
    return grid.size();
}
