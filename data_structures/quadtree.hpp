// quadtree.hpp

#include <array>
#include <vector>

#pragma once

namespace spatial {
    
    struct Range { long int start, end; };
    struct Rectangle { double xmin, xmax, ymin, ymax; };
    struct Point { double x, y; };

    Point midpoint(Rectangle rect) {
        double xmedian = (rect.xmin + rect.xmax) / 2;
        double ymedian = (rect.ymin + rect.ymax) / 2;
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
                    Node* parent;
                    int depth;
                    long long int code;
                    Range leaf_range;
                    Rectangle bounds;
                    Point center;
                    std::array<Node*,4> children; // [ *NW, *NE, *SW, *SE ]
                    Node(Node* parent, int depth, long long int code, Rectangle bounds);
                    void create_children();
            };

            struct Leaf {
                Node* node;
                std::vector<Datum<T>> bucket;
            };

            Node* root;
            std::vector<Leaf> leaves;

            Range recursive_build(Node* node, std::vector<Datum<T>> data);
            int get_quadrant(Point origin, Point p);

        public:
            Quadtree(double x0, double x1, double y0, double y1);
            Range build(std::vector<T>);
            int num_leaves();
    }; 
}
