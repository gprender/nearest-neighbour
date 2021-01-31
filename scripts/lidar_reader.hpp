// lidar_reader.hpp

#include <array>
#include <vector>
#include <fstream>

std::vector<std::string> split_string(std::string s, std::string delim);

class LidarReader {
    private:
        std::vector<std::vector<double>> point_data;
        std::array<double, 3> min;
        std::array<double, 3> max;

    public:
        LidarReader(std::string filename);
        std::vector<std::vector<double>> get_point_data();
        std::array<double, 3> get_min();
        std::array<double, 3> get_max();
};