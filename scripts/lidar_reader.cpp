// lidar_reader.cpp
/**
 * Full disclosure, this is NOT how one should be parsing lidar files
 * for practical applications; this just reads the .txt output from
 * LAStools' las2txt. Ideally, you should read directly from the .las
 * file with some existing library, or write your own .las parser.
 */

#include "lidar_reader.hpp"

LidarReader::LidarReader(std::string filename) {
    std::ifstream infile (filename);
    std::string line;

    if (!infile) {
        std::cout << "Unable to open data file: \"" << filename << "\"\n";
        std::cout << "Exiting...\n";
        exit(1);
    }

    // Parse the file header lines, denoted by '%' at position 0
    while (getline(infile, line)) {
        if (line[0] != '%') break;

        if (line.substr(2, 9).compare("min x y z") == 0) {
            std::string coords = line.substr(11, line.size() - 11);
            std::vector<std::string> split_coords = split_string(coords, " ");
            min[0] = stod(split_coords[0]);
            min[1] = stod(split_coords[1]); 
            min[2] = stod(split_coords[2]);      
        } 
        else if (line.substr(2, 9).compare("max x y z") == 0) {
            std::string coords = line.substr(11, line.size() - 11);
            std::vector<std::string> split_coords = split_string(coords, " ");
            max[0] = stod(split_coords[0]);
            max[1] = stod(split_coords[1]);
            max[2] = stod(split_coords[2]); 
        }
    }

    // Parse the actual point data (format: "x y z")
    do { 
        std::vector<std::string> str_coords = split_string(line, " ");
        std::vector<double> coords;
        for(auto const& s : str_coords) {
            coords.push_back(stod(s));
        }
        point_data.push_back(coords);
    } while (getline(infile, line)); 
}

std::vector<std::vector<double>> LidarReader::get_point_data() { return point_data; }

std::array<double, 3> LidarReader::get_min() { return min; }

std::array<double, 3> LidarReader::get_max() { return max; }

// This function is bad
std::vector<std::string> split_string(std::string s, std::string delim) {
    std::vector<std::string> v;

    int begin = s.find_first_not_of(delim);
    std::string rem = s.substr(begin, s.size() - begin);

    int next_split;
    std::string token;
    for (int i=0; i<3; i++) {
        next_split = rem.find_first_of(delim);
        token = rem.substr(0, next_split);
        rem = rem.substr(next_split + 1, rem.size() - next_split);
        v.push_back(token);
    }

    return v;
}