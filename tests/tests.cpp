// tests.cpp

#include <iostream>
#include <iomanip>
#include <algorithm>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "../data_structures/quadtree.cpp"
#include "lidar_reader.cpp"

/**
 * reg2048.txt is artificially generated in a 16x16 grid,
 * where each 1x1 cell contains 8 uniformly distributed points.
 * For a leaf size threshold of 8, we would expect this data 
 * to result in a complete quadtree of depth 4.
 */
std::string const reg2048 = "data/reg2048.txt";

/**
 * rand100k.txt contains 100,000 uniformly distributed points
 * inside of a square axis-aligned bounding box. This data is more
 * realistic than the above file, while still resulting in a
 * reasonably balanced quadtree.
 */
std::string const rand100k = "data/rand100k.txt";

using coord_t = spatial::coord_t;

/**
 * Verify that the k-nearest neighbours are ordered far -> close
 * This might seem like a weird property to aim for, but it helps verify
 * that the point priority queue used in queries is working as expected
 */
bool check_ordering(
    std::vector<std::vector<coord_t>> knn, 
    spatial::Point query_point
) {
    for (unsigned i=0; i<knn.size()-1; i++) {
        coord_t const dista = spatial::distance(
            query_point, (spatial::Point){knn[i][0], knn[i][1]}
        );
        coord_t const distb = spatial::distance(
            query_point, (spatial::Point){knn[i+1][0], knn[i+1][1]}
        );
        if (dista < distb) return false;
    }
    return true;
}

/**
 * Verify that we've correctly identified the k-nearest neighbours 
 * We do this with the following brute-force approach:
 * 
 *  1. Identify the point in the k-nn query result which is furthest 
 *      from the original query point (easy w/ the ordering property)
 * 
 *  2. Perform a linear scan through the entire point cloud
 * 
 *  3. If we find a point which is at least as close to the query point as the one
 *      we identified in (1.), we verify that point is a member of the query result
 * 
 * There's some ambiguity here if we come across the case where multiple points
 * in the quadtree are equidistant to the k'th nearest neighbour, in that we 
 * don't run any checks on such points. This might be bad, but it shouldn't
 * give false positives or false negatives for correctness. Be aware though!
 */
bool check_knn(
    std::vector<std::vector<coord_t>> knn, 
    spatial::Point query_point,
    std::vector<std::vector<coord_t>> point_data
) {
    coord_t const max_knn_dist = spatial::distance(
        query_point, (spatial::Point){knn[0][0], knn[0][1]}
    );
    for (auto const& point : point_data) {
        coord_t const current_dist = spatial::distance(
            query_point, (spatial::Point){point[0], point[1]}
        );
        // NOTE: We're not checking when current == max
        if (current_dist < max_knn_dist) {
            if (std::find(
                    point_data.begin(), 
                    point_data.end(), 
                    point
                ) == point_data.end()
            ) { return false; }
        }
    }
    return true;
}

TEST_CASE("Make sure all this quadtree code actually works", "[quadtree]") {

    SECTION("construction (point partitioning & recursion)") {
        LidarReader reader = LidarReader(reg2048);
        auto min = reader.get_min();
        auto max = reader.get_max();
        auto qt = new spatial::Quadtree<std::vector<coord_t>>(
            min[0], max[0], min[1], max[1]
        );
        qt->build(reader.get_point_data());

        REQUIRE(qt->num_leaves() == (16*16));
        REQUIRE(qt->check_depth_equals(4));

        delete qt;
    }

    SECTION("k-nearest neighbour querying") {
        LidarReader reader = LidarReader(rand100k);
        auto min = reader.get_min();
        auto max = reader.get_max();
        auto qt = new spatial::Quadtree<std::vector<coord_t>>(
            min[0], max[0], min[1], max[1]
        );
        auto point_data = reader.get_point_data();
        qt->build(point_data);

        auto knn1 = qt->query_knn(1, 100, 150);
        auto knn16 = qt->query_knn(16, 300, 450);
        auto knn32 = qt->query_knn(32, 250, 250);
        auto knnNW = qt->query_knn(8, 0, 0);
        auto knnSE = qt->query_knn(8, 500, 500);
        auto knnXX = qt->query_knn(16, 250, 750);

        REQUIRE(check_ordering(knn1, {100, 150}));
        REQUIRE(check_ordering(knn16, {300, 450}));
        REQUIRE(check_ordering(knn32, {250, 250}));
        REQUIRE(check_ordering(knnNW, {0, 0}));
        REQUIRE(check_ordering(knnSE, {500, 500}));
        REQUIRE(check_ordering(knnXX, {250, 750}));

        REQUIRE(check_knn(knn1, {100, 150}, point_data));
        REQUIRE(check_knn(knn16, {300, 450}, point_data));
        REQUIRE(check_knn(knn32, {250, 250}, point_data));
        REQUIRE(check_knn(knnNW, {0, 0}, point_data));
        REQUIRE(check_knn(knnSE, {500, 500}, point_data));
        REQUIRE(check_knn(knnXX, {250, 750}, point_data));

        delete qt;
    }
}
