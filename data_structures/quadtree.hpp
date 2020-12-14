// quadtree.hpp

#include <array>
#include <vector>
#include <memory>
#include <queue>

#include "spatial.hpp"

#pragma once

namespace spatial {

    template<typename T>
    class Quadtree {
        private:
            class Node {
                public:
                    int depth;
                    code_t code;
                    Rectangle bounds;
                    Point center;
                    Range leaf_range;
                    std::array<std::unique_ptr<Node>,4> children;

                    Node(int depth, code_t code, Rectangle bounds);
                    ~Node();
                    Range insert(
                        Quadtree<T>& tree, 
                        std::vector<Datum<T>> data
                    );
                    int get_quadrant(Point const p) const;
                    void create_children();
                    bool is_leaf() const;
                    
            };

            class NodePQ {
                private:
                    struct Element {
                        Node* node;
                        coord_t dist;
                    };

                    struct Closer {
                        bool operator()(Element const a, Element const b) {
                            return (a.dist > b.dist);
                        }
                    };

                    std::priority_queue<
                        Element, 
                        std::vector<Element>, 
                        Closer
                    > pq;
                    Point origin;

                public:
                    NodePQ(Point p):
                        origin(p)
                    { }

                    void push(Node* n) {
                        pq.push((Element){n, distance(origin, n->bounds)});
                    }

                    Element pop() {
                        auto const top_element = pq.top();
                        pq.pop();
                        return top_element;
                    }

                    Element const& peek() const { return pq.top(); }

                    // Assumes n->children is a collection of smart pointers
                    void expand(Node* n) {
                        for (auto const& child_ptr : n->children) {
                            push(child_ptr.get());
                        }
                    }
            };

            class DatumPQ {
                private:
                    struct Element {
                        Datum<T> datum;
                        coord_t dist;
                    };

                    struct Farther {
                        bool operator()(Element const a, Element const b) {
                            return (a.dist < b.dist);
                        }
                    };

                    std::priority_queue<
                        Element,
                        std::vector<Element>, 
                        Farther
                    > pq;
                    Point origin;

                public:
                    DatumPQ(Point p):
                        origin(p)
                    { }

                    void push(Datum<T> const& d) {
                        pq.push((Element){d, distance(origin, d.point)});
                    }

                    Element pop() {
                        auto const top_element = pq.top();
                        pq.pop();
                        return top_element;
                    }

                    Element const& peek() const { return pq.top(); }

                    void choose(Datum<T> const& d) {
                        coord_t const new_dist = distance(origin, d.point);
                        if (peek().dist > new_dist) {
                            pq.pop();
                            pq.push((Element){d, new_dist});
                        }
                    }

                    unsigned size() { return pq.size(); }

                    bool empty() { return (pq.empty()); }
            };

            struct Leaf { std::vector<Datum<T>> bucket; };

            // Node Priority Queue Element
            struct NodePQE {
                std::unique_ptr<Node> node;
                coord_t dist;
            };

            // Datum Priority Queue Element
            // Datum can be space-intensive, should use a pointer here instead
            struct DatumPQE {
                Datum<T> datum;
                coord_t dist;
            };
            
            std::unique_ptr<Node> root;
            std::vector<std::vector<Datum<T>>> leaves;

        public:
            Quadtree(coord_t x0, coord_t x1, coord_t y0, coord_t y1);
            ~Quadtree();
            Range build(std::vector<T> const& raw_data);
            std::vector<T> query_knn(
                unsigned const k, coord_t const x, coord_t const y
            ) const;
            int num_leaves() const;
            bool depth_equals(unsigned const k) const;
    }; 
}
