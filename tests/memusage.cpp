// memusage.cpp
/**
 * Usage: ./memusage <data1.txt> <data2.txt> ...
 * 
 * Note that this testing program uses windows-specific header files,
 * and as such will not work on Linux machines. Requires -lpsapi when compiling
 */

#include <iostream>
#include <iomanip>

#include "windows.h"
#include "psapi.h"

#include "../data_structures/quadtree.cpp"
#include "../data_structures/rtree.cpp"
#include "lidar_reader.cpp"

using coord_t = spatial::coord_t;

size_t process_memusage() {
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(
        GetCurrentProcess(), 
        (PROCESS_MEMORY_COUNTERS*)&pmc, 
        sizeof(pmc)
    );
    return pmc.WorkingSetSize;
} 

void quadtree_benchmark(std::string const filename) {
    std::cout << "\nRunning quadtree memory benchmark for \'" 
              << filename << "\':\n";

    // Parsing the data & header
    LidarReader reader(filename);
    auto const& min = reader.get_min();
    auto const& max = reader.get_max();

    // Record a memory baseline here
    // The space usage of reading the data isn't relevant 
    auto mem_baseline = process_memusage();

    spatial::Quadtree<std::vector<coord_t>> qt(
        min[0], max[0], min[1], max[1]
    );

    std::cout << "\tBuilding the quadtree... \t";
    qt.build(reader.get_point_data());
    std::cout << process_memusage() - mem_baseline << " bytes used\n";
}

void rtree_benchmark(std::string filename) {
    std::cout << "\nRunning R-tree memory benchmark for \'" 
              << filename << "\':\n";

    LidarReader reader = LidarReader(filename);

    auto mem_baseline = process_memusage();

    std::cout << "\tBuilding the R-tree...   \t";
    spatial::Rtree<std::vector<double>> rtree;
    rtree.build(reader.get_point_data());

    std::cout << process_memusage() - mem_baseline << " bytes used\n";
}

int main(int argc, char** argv) {

    for (int i=1; i<argc; i++) {
        quadtree_benchmark(argv[i]);
        rtree_benchmark(argv[i]);
    }

    return 0;
}
