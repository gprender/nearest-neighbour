// quadtree.hpp

#include <array>
#include <vector>

#pragma once

// Type aliases! These could also go inside of the namespace
using coord_t = double; // Needs to fit whatever numbers are used for data
using code_t = long long int; // Needs at least 2*(quadtree height) bits
using index_t = unsigned long int; // Needs to fit the total # of leaves

namespace spatial {
 
    struct Range { index_t start, end; };
    struct Rectangle { coord_t xmin, xmax, ymin, ymax; };
    struct Point { coord_t x, y; };

    Point midpoint(Rectangle const rect) {
        coord_t const xmedian = (rect.xmin + rect.xmax) / 2;
        coord_t const ymedian = (rect.ymin + rect.ymax) / 2;
        return {xmedian, ymedian};
    }

    /**
     * A "datum" is a single element in the quadtree. Contains the raw data T,
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
            
            Node* root;
            std::vector<Leaf> leaves;

            Range recursive_build(Node* const node, std::vector<Datum<T>> data);
            void recursive_deconstruct(Node* const node);
            int get_quadrant(Point const origin, Point const p);
            Node* find_leaf(Point const p);
            Node* traverse(Node* const origin_node, code_t const target_code);
            std::vector<code_t> calc_neighbour_codes(
                code_t const origin_code, int const depth);

        public:
            Quadtree(coord_t x0, coord_t x1, coord_t y0, coord_t y1);
            ~Quadtree();
            Range build(std::vector<T> const raw_data);
            std::vector<T> query(coord_t const x, coord_t const y);
            int num_leaves();

            /**
             * A lightweight finite state machine for computing neighbour codes
             * First dimension is a 2-bit code digit (0-3)
             * Second dimension is the direction (N, NE, E, SE, S, SW, W, NW)
             * The third dimension is a <digit, halt?> pair: "digit" is the new
             * code digit after traveling in the given direction, and "halt?"
             * is either a signal to stop the calculation because we're at a
             * sibling (-1), or the new travel direction (0-7).
             */
            int const fsm[4][8][2] = {
                {{2,0}, {3,0}, {1,-1}, {3,-1}, {2,-1}, {3,6}, {1,6}, {3,7}},
                {{3,0}, {2,1}, {0,2}, {2,2}, {3,-1}, {2,-1}, {0,-1}, {2,0}},
                {{0,-1}, {1,-1}, {3,-1}, {1,4}, {0,4}, {1,5}, {3,6}, {1,6}},
                {{1,-1}, {0,2}, {2,2}, {0,3}, {1,4}, {0,4}, {2,-1}, {0,-1}}
            };
    }; 
}
