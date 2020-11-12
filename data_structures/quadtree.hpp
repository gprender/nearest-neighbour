// quadtree.hpp

#include <array>
#include <vector>

#include "spatial.hpp"

#pragma once

namespace spatial {

    template<typename T>
    class Quadtree {
        private:
            class Node {
                public:
                    int depth;
                    Node* parent;
                    code_t code;
                    Rectangle bounds;
                    Point center;
                    Range leaf_range;
                    std::array<Node*,4> children; // [ *NW, *NE, *SW, *SE ]

                    Node(Node* parent, int depth, code_t code, Rectangle bounds);
                    void create_children();
                    bool is_leaf();
            };

            struct Leaf {
                Node* node;
                std::vector<Datum<T>> bucket;
            };

            // Node Priority Queue Element
            struct NodePQE {
                Node* node;
                coord_t dist;
            };

            // Datum Priority Queue Element
            // Datum can be space-intensive, should use a pointer here instead
            struct DatumPQE {
                Datum<T> datum;
                coord_t dist;
            };
            
            Node* root;
            std::vector<Leaf> leaves;

            Range recursive_build(
                Node* const node, std::vector<Datum<T>> data
            );
            void recursive_deconstruct(Node* const node);
            int get_quadrant(Point const origin, Point const p) const;

        public:
            Quadtree(coord_t x0, coord_t x1, coord_t y0, coord_t y1);
            ~Quadtree();
            Range build(std::vector<T> const raw_data);
            std::vector<T> query_knn(
                unsigned const k, coord_t const x, coord_t const y
            ) const;
            int num_leaves() const;
            bool depth_equals(unsigned const k) const;
    }; 
}
