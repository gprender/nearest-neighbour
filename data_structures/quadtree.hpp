// quadtree.hpp

#include <vector>

#pragma once

namespace spatial {

    struct Point { double x, y; };
    struct Rectangle { double xmin, xmax, ymin, ymax; };

    /**
     * A "datum" is a single element in the quadtree. Contains the raw data
     * passed in by a call to insert(), as well as a 2d interpretation of that
     * data as a "point".
     */
    template<typename T>
    struct Datum {
        T data;
        Point point;
    };

    template<typename T>
    class Quadtree {
        private:
            int height, load;
            struct Rectangle bounds;
            double xmedian, ymedian;
            Datum<T>* data; // Collection of elements in a leaf - is there a better name for this?
            Quadtree *NW, *NE, *SW, *SE;
            int split_node();
            int recursive_insert(Datum<T>);

        public:
            Quadtree(double x0, double x1, double y0, double y1);
            int insert(T);
            int insert(Datum<T>);
            int getHeight();
            int getLoad();
            bool isLeaf();
    }; 

}