// timing.cpp

#include <iostream>
#include <iomanip>
#include <chrono>

#include "../data_structures/quadtree.cpp"
#include "lidar_reader.cpp"

std::string const DF1 = "data/rand100k.txt";
std::string const DF2 = "data/rand500k.txt";
std::string const DF3 = "data/rand1m.txt";
std::string const QDF = "data/rand1k.txt";

using coord_t = spatial::coord_t;

void benchmark(std::string filename) {
    std::cout << "\nRunning timing benchmark for \'" << filename << "\':\n";

    // Parsing the data & header
    LidarReader reader = LidarReader(filename);
    auto const min = reader.get_min();
    auto const max = reader.get_max();

    auto qt = new spatial::Quadtree<std::vector<coord_t>>(
        min[0], max[0], min[1], max[1]
    );

    // Quadtree construction
    std::cout << "\tBuilding the quadtree... ";
    auto start = std::chrono::system_clock::now();
    qt->build(reader.get_point_data());
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>
        (end - start).count();
    std::cout << elapsed << " microseconds\n";

    // k-nn queries
    LidarReader qreader = LidarReader(QDF);
    std::cout << "\tQuerying k-nearest neighbours x1000...\n";
    for (auto k : {1, 8, 16, 32}) {
        std::cout << "\t\tk=" << k << ":\t";
        start = std::chrono::system_clock::now();
        for (auto e : qreader.get_point_data()) {
            qt->query_knn(k, e[0], e[1]);
        }
        end = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::microseconds>
            (end - start).count();
        std::cout << elapsed << " microseconds\n";
    }

    std::cout << "\n";

    delete qt;
}

int main() {

    benchmark(DF1);

    benchmark(DF2);

    benchmark(DF3);

    return 0;
}
