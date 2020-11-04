// quadtree.hpp

#include <cmath>
#include <array>
#include <vector>

#pragma once

namespace spatial {
    // Type aliases!
    using coord_t = double; // Needs to fit whatever numbers are used for data
    using code_t = long long int; // Needs at least 2*(quadtree height) bits
    using index_t = unsigned long int; // Needs to fit the total # of leaves
 
    struct Range { index_t start, end; };
    struct Rectangle { coord_t xmin, xmax, ymin, ymax; };
    struct Point { coord_t x, y; };

    Point midpoint(Rectangle const rect) {
        coord_t const xmedian = (rect.xmin + rect.xmax) / 2;
        coord_t const ymedian = (rect.ymin + rect.ymax) / 2;
        return {xmedian, ymedian};
    }

    coord_t distance(Point const p, Point const q) {
        coord_t const dx = p.x - q.x;
        coord_t const dy = p.y - q.y;
        return std::sqrt(dx*dx + dy*dy);
    }

    coord_t distance(Point const p, Rectangle const rect) {
        coord_t dx = std::max(rect.xmin - p.x, p.x - rect.xmax);
        coord_t dy = std::max(rect.ymin - p.y, p.y - rect.ymax);
        dx = std::max(dx, 0.0);
        dy = std::max(dy, 0.0);
        return std::sqrt(dx*dx + dy*dy);
    }

    /**
     * A "datum" is a single element in the quadtree. Contains the raw data,
     * as well as a 2d interpretation of that data as a "point".
     */
    template<typename T>
    struct Datum {
        T data;
        Point point;
    };

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
