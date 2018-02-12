#ifndef CLASSES_HPP_INCLUDED
#define CLASSES_HPP_INCLUDED

#include <vector>
#include <string.h>

using namespace std;

const float PI = 3.141592653589793238463;

typedef std::vector< std::vector<double> > vvec;
typedef vector< vector<int> > vint;

/** command line arguments **/

struct CommandLine {

    string file_path;                 // -i option, input file
    string output_path;               // -o option, output file
    string lower_slice;               // -l option, lower height
    string upper_slice;               // -u option, upper height
    double  pixel_size;               // -p option, pixel size in meters
    double  max_radius;               // -r option, maximum radius to test
    double  min_density;              // -d option, minimum density to consider on the Hough transform
    int  min_votes;                   // -v option, minimum votes count at the output
    string output_las;                // -O option, save a las/laz/txt output
    bool help;                        // -h/? option, get help
    float height_interval;            // -z option, z interval to slice tree
    string output_stack;
    string output_stack_coordinates;
    string clouds_directory;
    string reports_directory;

};

class CloudStats{
    public:
        double x_min;
        double x_max;
        double y_min;
        double y_max;
        double z_min;
        double z_max;
        int point_count;
        double points_per_zm;

        void getPPM(){
            points_per_zm = point_count / (z_max - z_min);
        }
};

class Slice{
    public:
        vvec slice;
        vector<double> dims;
};

class Raster{

    public:
        vint matrix;
        int x_dim;
        int y_dim;
        int max_count;
        float pixel_size;
        double thickness;
        double min_x;
        double max_x;
        double min_y;
        double max_y;
        double min_z;
        double max_z;

        void setMatrixSize(int nx, int ny){
            matrix.resize(nx, vector<int>(ny));
        };

};

class HoughCircle{
    public:
        float x_center;
        float y_center;
        double radius;
        int n_votes;
};


class HoughCenters{
    public:
        vector<HoughCircle> circles;
        int main_circle;
        double avg_x;
        double avg_y;
        float aggregate_radius;
        double low_z;
        double up_z;

        void getCenters(){
            double ax = 0, ay = 0;
            main_circle = 0;
            for(unsigned i = 0; i < circles.size(); ++i){
                ax += circles[i].x_center;
                ay += circles[i].y_center;

                if(circles[i].n_votes > circles[main_circle].n_votes){
                    main_circle = i;
                }
            }
            avg_x = ax / circles.size();
            avg_y = ay / circles.size();
        }

        HoughCircle getMainCircle(){
            return circles[main_circle];
        }
};

class StemSegment{
    public:
       HoughCircle model_circle;
       double z_min;
       double z_max;
       int n_points;
};

class PlotCloudTreeSlice{
    public:
        unsigned int n_tree;
        unsigned int n_slice;
        float x_average;
        float y_average;
        float x_main;
        float y_main;
        float radius;
        float z_min;
        float z_max;
        unsigned int votes;
};

#endif // CLASSES_HPP_INCLUDED
