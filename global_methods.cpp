#include <iostream>
#include "lasreader.hpp"
#include "laswriter.hpp"
#include "lasfilter.hpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <unordered_map>
#include <set>
#include <vector>
#include <typeinfo>
#include <fstream>
#include <algorithm>
#include <string.h>
#include <sstream>
#include "classes.hpp"
#include "plot_methods.hpp"

using namespace std;

//output suffix
string outputNameAppend(string path, string suffix){

    unsigned last_dot = path.find_last_of(".");
    string new_path = path.substr(0, last_dot) + suffix;

    return new_path;
}

//get cloud statistics
CloudStats getStats(string file){

  LASreadOpener lasreadopener;
  lasreadopener.set_file_name(file.c_str());
  LASreader* lasreader = lasreadopener.open();

    unsigned last_dot = file.find_last_of(".");
    string file_format = file.substr(last_dot+1, file.length());

    CloudStats result;

    if(file_format == "las" || file_format == "laz"){

    result.x_min = lasreader->get_min_x();
    result.x_max = lasreader->get_max_x();
    result.y_min = lasreader->get_min_y();
    result.y_max = lasreader->get_max_y();
    result.z_min = lasreader->get_min_z();
    result.z_max = lasreader->get_max_z();
    result.point_count = lasreader->npoints;
    result.getPPM();

    } else {

    unsigned i = 1;
    double x, y, z, min_x, max_x, min_y, max_y, min_z, max_z;

    while(lasreader->read_point()){

        x = lasreader->get_x();
        y = lasreader->get_y();
        z = lasreader->get_z();

        if(i == 1){
           min_x = max_x = x;
           min_y = max_y = y;
           min_z = max_z = z;
        }

        if(x > max_x) max_x = x;
        if(y > max_y) max_y = y;
        if(x < min_x) min_x = x;
        if(y < min_y) min_y = y;
        if(z < min_z) min_z = z;
        if(z > max_z) max_z = z;

        i++;

    }

    result.x_min = min_x;
    result.x_max = max_x;
    result.y_min = min_y;
    result.y_max = max_y;
    result.z_min = min_z;
    result.z_max = max_z;
    result.point_count = i;
    result.getPPM();

    }

    lasreader->close();
    delete lasreader;

    cout <<
        "# x: " << result.x_min << " : " << result.x_max << "\n" <<
        "# y: " << result.y_min << " : " << result.y_max << "\n" <<
        "# z: " << result.z_min << " : " << result.z_max << "\n" <<
        "# points: " << result.point_count << " : " << result.points_per_zm << " /m above ground"
    << endl;

    return result;

};

//calculate absolute center coordinate based on pixel position
vector<double> absCenter(int x, int y, double min_x, double min_y, float step){
    double x_cen = ( min_x + (step/2) ) + ( x * step );
    double y_cen = ( min_y + (step/2) ) + ( y * step );

    return vector<double> {x_cen, y_cen};
};

//calculate pixel based on absolute coordinate
vector<int> pixPosition(double x, double y, double min_x, double min_y, float step){
    int x_pix = floor( (x - min_x) / step );
    int y_pix = floor( (y - min_y) / step );

    return vector<int> {x_pix, y_pix};
};
