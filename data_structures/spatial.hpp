// spatial.hpp
/**
 * A header file containing some simple 2d structures and functions
 */

#include <cmath>

#pragma once

namespace spatial {

    // Type aliases!
    using coord_t = double;  // Needs to fit whatever numbers are used for data
    using area_t = coord_t;  // Better semantics for area calculations
    using code_t = long long int;  // Needs at least 2*(quadtree height) bits
    using index_t = unsigned long int;  // Needs to fit the total # of leaves
 
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

    area_t area(Rectangle const rect) {
        return (rect.xmax - rect.xmin)*(rect.ymax - rect.ymin);
    }

    /**
     * Calculate the minimum bounding box of a rectangle and a point.
     */
    Rectangle min_bounding_box(Rectangle const rect, Point const p) {
        return (Rectangle) {
            std::min(rect.xmin, p.x),
            std::max(rect.xmax, p.x),
            std::min(rect.ymin, p.y),
            std::max(rect.ymax, p.y)
        };
    }

    /**
     * Calculate the minimum bounding box of two rectangles.
     */
    Rectangle min_bounding_box(Rectangle const r1, Rectangle const r2) {
        return (Rectangle) {
            std::min(r1.xmin, r2.xmin),
            std::max(r1.xmax, r2.xmax),
            std::min(r1.ymin, r2.ymin),
            std::max(r1.ymax, r2.ymax)
        };
    }

    /**
     * Verify whether one rectangle contains another.
     */
    bool contains(Rectangle const outer, Rectangle const inner) {
        return (
            outer.xmin <= inner.xmin 
            && outer.xmax >= inner.xmax 
            && outer.ymin <= inner.ymin 
            && outer.ymax >= inner.ymax
        );
    }

    bool contains(Rectangle const rect, Point const p) {
        return (
            rect.xmin <= p.x
            && rect.xmax >= p.x 
            && rect.ymin <= p.y 
            && rect.ymax >= p.y
        );
    }

    void print_rect(Rectangle const rect) {
        std::cout << std::fixed << std::setprecision(2) 
                  << "x[" << rect.xmin << ", " << rect.xmax << "]  ";

        std::cout << std::fixed << std::setprecision(2) 
                  << "y[" << rect.ymin << ", " << rect.ymax << "]\n";
    }

    /**
     * A single element in a 2d space partitioning tree.
     * Contains raw data, and an interpetation of that data as a 2d point.
     */
    template<typename T>
    struct Datum {
        T data;
        Point point;
    };

}
