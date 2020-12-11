// timing.cpp
/**
 * Usage: ./timing <data1.txt> <query1.txt> <data2.txt> <query2.txt> ...
 * Where 'data_.txt' contains the point data for the tree constructions,
 * and 'query_.txt' contains the points which we will query around.
 */

#include <iostream>
#include <iomanip>
#include <chrono>

#include "../data_structures/quadtree.cpp"
#include "../data_structures/rtree.cpp"
#include "lidar_reader.cpp"

using coord_t = spatial::coord_t;

void quadtree_benchmark(std::string data_file, std::string query_file) {
    std::cout << "\nRunning quadtree timing benchmark for \'" << data_file << "\',\n"
              << "using " << query_file << " for query points.\n";

    // Parsing the data & header
    LidarReader reader(data_file);
    auto const& min = reader.get_min();
    auto const& max = reader.get_max();

    spatial::Quadtree<std::vector<coord_t>> qt(
        min[0], max[0], min[1], max[1]
    );

    // Quadtree construction
    std::cout << "\tBuilding the quadtree... ";
    auto start = std::chrono::system_clock::now();
    qt.build(reader.get_point_data());
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
        (end - start).count();
    std::cout << elapsed << " milliseconds\n";

    // k-nn queries
    LidarReader query_reader(query_file);
    std::cout << "\tQuerying k-nearest neighbours x1000...\n";
    for (auto const k : {1, 8, 32}) {
        std::cout << "\t\tk=" << k << ":\t";
        coord_t filler = 0;
        start = std::chrono::system_clock::now();
        for (auto const& p : query_reader.get_point_data()) {
            auto const knn = qt.query_knn(k, p[0], p[1]);
            filler += knn[0][2];
        }
        end = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
            (end - start).count();
        std::cout << elapsed << " milliseconds";
        std::cout << "  \t(filler: " << filler << ")\n";
    }
    std::cout << "\n";
}

void rtree_benchmark(std::string data_file, std::string query_file) {
    std::cout << "\nRunning R-tree timing benchmark for \'" << data_file << "\',\n"
              << "using " << query_file << " for query points.\n";

    LidarReader data_reader(data_file);
    spatial::Rtree<std::vector<double>> rtree;

    std::cout << "\tBuilding the R-tree... ";
    auto start = std::chrono::system_clock::now();
    rtree.build(data_reader.get_point_data());
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
        (end - start).count();
    std::cout << elapsed << " milliseconds\n";

    LidarReader query_reader(query_file);
    std::cout << "\tQuerying k-nearest neighbours x1000...\n";
    for (auto const k : {1, 8, 32}) {
        std::cout << "\t\tk=" << k << ":\t";
        coord_t filler = 0;
        start = std::chrono::system_clock::now();
        for (auto const& p : query_reader.get_point_data()) {
            auto const knn = rtree.query_knn(k, p[0], p[1]);
            filler += knn[0][2];
        }
        end = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
            (end - start).count();
        std::cout << elapsed << " milliseconds";
        std::cout << "  \t(filler: " << filler << ")\n";
    }
    std::cout << "\n";
}

int main(int argc, char** argv) {

    for (int i=1; i<argc; i+=2) {
        quadtree_benchmark(argv[i], argv[i+1]);
        rtree_benchmark(argv[i], argv[i+1]);
    }

    return 0;
}
