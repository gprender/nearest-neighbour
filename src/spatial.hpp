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

    void print_range(Range const r) {
        std::cout << "[ " << r.start << ", " << r.end << " ]\n";
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

    /**
     * Build a vector of Datum<T> from a vector of <T>.
     */
    template<typename T>
    std::vector<Datum<T>> datumize(std::vector<T> const& raw_data) {
        std::vector<Datum<T>> formatted_data;
        formatted_data.reserve(raw_data.size());
        for (auto const& raw_datum : raw_data) {
            Point const new_point = {raw_datum[0], raw_datum[1]};
            Datum<T> const new_datum = {raw_datum, new_point};
            formatted_data.push_back(new_datum);
        }
        return formatted_data;
    }

    /**
     * For a 1d range [min,max] divided into 'dim' equal partitions, 
     * find the partition (or index) which contains 'coord'.
     */
    int grid_index(
        coord_t const coord, 
        coord_t const min, 
        coord_t const max, 
        int const dim
    ) {
        return (coord - min) * dim / (max - min);
    }

    // Bitstring constants for use in interleaving integers
    uint8_t const _shifts[] = { 1, 2, 4, 8 };
    uint32_t const _masks[] = { 
        0x55555555, 
        0x33333333, 
        0x0F0F0F0F, 
        0x00FF00FF 
    };

    /**
     * Space 16-bits out into 32-bits, with zeros inbetween.
     * e.g., 1111 -> 01010101
     */
    uint32_t space_bits(uint16_t const i0) {
        uint32_t i = (i0 | (i0 << _shifts[3])) & _masks[3];
        i = (i | (i << _shifts[2])) & _masks[2];
        i = (i | (i << _shifts[1])) & _masks[1];
        return (i | (i << _shifts[0])) & _masks[0];
    }
    
    /**
     * Interleave two 16-bit integers into a 32-bit integer.
     * e.g., (ABCD, EFGH) -> EAFB GCHD
     */
    uint32_t interleave(uint16_t const a, uint16_t const b) {
        return space_bits(a) | (space_bits(b) << 1);
    }

}
