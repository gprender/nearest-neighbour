// memusage.cpp
/**
 * Note that this testing program uses windows-specific header files,
 * and as such will not work on Linux machines. Requires -lpsapi when compiling
 */

#include <iostream>
#include <iomanip>

#include "windows.h"
#include "psapi.h"

#include "../data_structures/quadtree.cpp"
#include "lidar_reader.cpp"

std::string const DF1 = "data/rand100k.txt";
std::string const DF2 = "data/rand500k.txt";
std::string const DF3 = "data/rand1m.txt";

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

void benchmark(std::string filename) {
    std::cout << "\nRunning memory benchmark for \'" << filename << "\':\n";

    // Parsing the data & header
    LidarReader reader = LidarReader(filename);
    auto const min = reader.get_min();
    auto const max = reader.get_max();

    // Record a memory baseline here
    // The space usage of reading the data isn't relevant 
    auto mem_baseline = process_memusage();

    auto qt = new spatial::Quadtree<std::vector<coord_t>>(
        min[0], max[0], min[1], max[1]
    );

    std::cout << "\tBuilding the quadtree... ";
    qt->build(reader.get_point_data());
    std::cout << process_memusage() - mem_baseline << " bytes used\n";

    delete qt;
}

int main() {

    benchmark(DF1);

    benchmark(DF2);

    benchmark(DF3);

    return 0;
}
