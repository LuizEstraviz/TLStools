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

/*
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
*/
using namespace std;

/*** types ***/
typedef vector< vector<double> > vvec;
typedef vector< vector<int> > vint;

const float PI = 3.141592653589793238463;

/*** classes ***/
class Slice{
    public:
        vvec slice;
        array<double,6> dims;
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

class Circle{
    public:
        vvec perimeter_points;
        vector< double > radii;
        float max_radius;
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
};
/*** functions ***/

/*
//string conversion
std::string to_string(int value){

      //create an output string stream
      std::ostringstream os ;

      //throw the value into the string stream
      os << value ;

      //convert the string stream into a string and return
      return os.str() ;

    }

//string split
vector<string> str_split(string str){

    string buf; // Have a buffer string
    stringstream ss(str); // Insert the string into a stream

    vector<string> tokens; // Create vector to hold our words

    while (ss >> buf){
        //int bint = stoi(buf);
        tokens.push_back(buf);
    }

    return tokens;
}
*/

//output suffix
string outputNameAppend(string path, string suffix = "_result.txt"){

    unsigned last_dot = path.find_last_of(".");
    string new_path = path.substr(0, last_dot) + suffix;

    return new_path;
}

//returns all points in between two heights and their x,y range
Slice getSlice(string file, string lower = "1.0", string upper = "2.0"){

    Slice slc;

    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(file.c_str());

    char* MY_other_argv[4];

    MY_other_argv[0] = (char*)malloc(1000); strcpy(MY_other_argv[0], " ");
    MY_other_argv[1] = (char*)malloc(1000); strcpy(MY_other_argv[1], "-keep_z");
    MY_other_argv[2] = (char*)malloc(1000); strcpy(MY_other_argv[2], lower.c_str() );
    MY_other_argv[3] = (char*)malloc(1000); strcpy(MY_other_argv[3], upper.c_str() );

    lasreadopener.parse(4, MY_other_argv);
    LASreader* lasreader = lasreadopener.open();

    int i = 0;
    double max_x;
    double max_y;
    double min_x;
    double min_y;
    double min_z;
    double max_z;

    //vvec xyz_matrix;
    vector<double> xyz_row(3);
    double x;
    double y;
    double z;

    while(lasreader->read_point()){

        x = lasreader->get_x();
        y = lasreader->get_y();
        z = lasreader->get_z();

        xyz_row[0] = x;
        xyz_row[1] = y;
        xyz_row[2] = z;

        slc.slice.push_back( xyz_row );

         if(i == 0){
           max_x = x;
           max_y = y;
           min_x = x;
           min_y = y;
           min_z = z;
           max_z = z;
         }

        if(x > max_x) max_x = x;
        if(y > max_y) max_y = y;
        if(x < min_x) min_x = x;
        if(y < min_y) min_y = y;
        if(z < min_z) min_z = z;
        if(z > max_z) max_z = z;

        i++;

    }

    lasreader->close();

    slc.dims = {min_x, max_x, min_y, max_y, min_z, max_z};

    return slc;

}

//returns a matrix of counts and maximum count value
Raster getCounts(Slice* slice , float pixel_size = 0.025){

    Raster ras;
    ras.min_x = (*slice).dims[0];
    ras.max_x = (*slice).dims[1];
    ras.min_y = (*slice).dims[2];
    ras.max_y = (*slice).dims[3];
    ras.min_z = (*slice).dims[4];
    ras.max_z = (*slice).dims[5];
    ras.pixel_size = pixel_size;
    ras.thickness = ras.max_z - ras.min_z;

    int xn = abs( ceil( (ras.max_x - ras.min_x) / pixel_size ) ) ;
    int yn = abs( ceil( (ras.max_y - ras.min_y) / pixel_size ) ) ;

    ras.x_dim = xn;
    ras.y_dim = yn;

    vint matrix( xn, vector<int>(yn) );
    int max_votes = 0;
    double point[3];
    int xCell;
    int yCell;

    for(unsigned int i = 0; i < (*slice).slice.size(); i++){

        point[0] = (*slice).slice[i][0];
        point[1] = (*slice).slice[i][1];
        point[2] = (*slice).slice[i][2];

        xCell = abs( floor( (point[0] - ras.min_x) / pixel_size ) );
        yCell = abs( floor( (point[1] - ras.min_y) / pixel_size ) );

        if(xCell >= xn) xCell = xn-1;
        if(yCell >= yn) yCell = yn-1;

        matrix[xCell][yCell] += 1;

        if(matrix[xCell][yCell] > max_votes) max_votes = matrix[xCell][yCell];

    }

    ras.matrix = matrix;
    ras.max_count = max_votes;

    return ras;

}
/*
Circle buildCircles(float r_min = 0.025, float r_max = 0.5){

    Circle circ;
    circ.max_radius = r_max;
    vector<double> radii;
    vector<float> angles;
    float angle_step;// = r_min / r_max;

    for(double i = r_min; i <= r_max; i += r_min){
        radii.push_back(i);
    }

    circ.radii = radii;

    vector<double> coord(3);
    for(unsigned i = 0; i < radii.size(); i++){

     angle_step = r_min / radii[i];
     for(float j = 0; j <= PI*2; j += angle_step){
        angles.push_back(j);
     }

        for(unsigned j = 0; j < angles.size(); j++){

            coord[0] = cos(angles[j])*radii[i];
            coord[1] = sin(angles[j])*radii[i];
            coord[2] = radii[i];

            circ.perimeter_points.push_back(coord);

        }
        angles = {};
    }

    return circ;
};
*/
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
/*
//sum center values to circle points
Circle sumCenters(Circle circle, double xc, double yc){
    Circle out_circle = circle;

    for(unsigned i=0; i < circle.perimeter_points.size(); i++){
        out_circle.perimeter_points[i][0] += xc;
        out_circle.perimeter_points[i][1] += yc;
    }

    return out_circle;
};
*/
vint rasterCircle(float radius, float pixel_size = 0.025, double cx = 0, double cy = 0, double mx = 0, double my = 0){

    int n_points = ceil( (2 * PI * radius) / pixel_size );
    double angle_dist = 2 * PI / n_points;
    double x, y;
    vector<int> pxy(2);
    vint pixels;

    for(double i = 0; i < 2*PI; i += angle_dist){
        x = cos(i)*radius + cx;
        y = sin(i)*radius + cy;

        pxy = pixPosition(x, y, mx, my, pixel_size);
        pixels.push_back(pxy);
    }

    return pixels;

};

vector<HoughCenters> getCenters(Raster* raster, float max_radius = 0.5, float min_den = 0.1, int min_votes = 2){

    //count raster properties (&raster)
    unsigned int x_len = raster->x_dim;
    unsigned int y_len = raster->y_dim;
    int min_count = floor( (raster->max_count)*min_den );

    //empty raster properties (votes)
    Raster empty_raster;
    empty_raster.min_x = raster->min_x - max_radius;
    empty_raster.min_y = raster->min_y - max_radius;
    empty_raster.max_x = raster->max_x + max_radius;
    empty_raster.max_y = raster->max_y + max_radius;
    empty_raster.pixel_size = raster->pixel_size;
    empty_raster.x_dim = abs( ceil( (empty_raster.max_x - empty_raster.min_x) / empty_raster.pixel_size ) );
    empty_raster.y_dim = abs( ceil( (empty_raster.max_y - empty_raster.min_y) / empty_raster.pixel_size ) );
    empty_raster.setMatrixSize(empty_raster.x_dim, empty_raster.y_dim);

    //get valid pixels
    vector< vector<unsigned int> > pixels;

    for(unsigned i = 0; i < x_len; ++i){
        for(unsigned j = 0; j < y_len; ++j){
            if(raster->matrix[i][j] >= min_count ){
                pixels.push_back( {i,j} );
            }
        }
    }

    //make circles centered in every valid pixel
    vector<HoughCenters> g_circles = {};
    HoughCenters g_cen;
    g_cen.circles = {};
    g_cen.avg_x = -1000;
    g_cen.avg_y = -1000;
    g_cen.low_z = raster->min_z;
    g_cen.up_z = raster->max_z;
    g_cen.aggregate_radius = max_radius * 3;
    g_circles.push_back(g_cen);

    vint votes;
    HoughCircle hc;
    vector<double> center;
    vint h_circle;
    set<unsigned long long int> pixel_set; //cod = 100.000*x + y
    vector<double> coor(2);
    //vector<HoughCircle> p_circles;
    unsigned int vx, vy;
    vector<int> pixel(2);
    vector<float> radii;

    for(float i = 0; i <= max_radius; i+=raster->pixel_size){
        radii.push_back(i);
    }

    for(auto& rad : radii){
    votes = empty_raster.matrix;
    hc.radius = rad;
    pixel_set = {};

        for(unsigned i = 0; i < pixels.size(); ++i){
            center = absCenter(pixels[i][0], pixels[i][1], raster->min_x, raster->min_y, raster->pixel_size);
            h_circle = rasterCircle(rad, empty_raster.pixel_size, center[0], center[1], empty_raster.min_x, empty_raster.min_y);

            for(unsigned j = 0; j < h_circle.size(); ++j){
                vx = h_circle[j][0];
                vy = h_circle[j][1];
                votes[ vx ][ vy ] += 1;

                if(votes[ vx ][ vy ] >= min_votes){
                    pixel_set.insert( 100000*vx + vy );
                }
            }
        }

        for(auto& k : pixel_set){
            vx = floor(k / 100000);
            vy = k - 100000*vx;
            coor = absCenter(vx, vy, empty_raster.min_x, empty_raster.min_y, empty_raster.pixel_size);
            hc.x_center = coor[0];
            hc.y_center = coor[1];
            hc.n_votes = votes[vx][vy];
            //p_circles.push_back(hc);

            for(unsigned i = 0; i < g_circles.size(); ++i){
                double dist = sqrt( pow(hc.x_center - g_circles[i].avg_x, 2) + pow(hc.y_center - g_circles[i].avg_y, 2) );

                if(dist < g_cen.aggregate_radius){
                    g_circles[i].circles.push_back(hc);
                    break;
                }

                if(i == (g_circles.size()-1)){
                    g_cen.circles = {};
                    g_cen.circles.push_back(hc);
                    g_cen.avg_x = hc.x_center;
                    g_cen.avg_y = hc.y_center;
                    g_circles.push_back(g_cen);
                }
            }
        }

    }
    return g_circles;

}

void getPreciseCenters(vector<HoughCenters>& circles){

    circles.erase(circles.begin());
    for(unsigned i = 0; i < circles.size(); ++i){

        circles[i].getCenters();

    }

};

void saveReport(vector<HoughCenters>& centers, string file_path = "result.txt"){

        ofstream result_file(file_path);

        result_file << "x_average   y_average   x_main   y_main   r_main   votes   z_min   z_max" << endl;

        for(unsigned int i = 0; i < centers.size() ; ++i){

            result_file << centers[i].avg_x << "    " << centers[i].avg_y << "   " <<
            centers[i].circles[ centers[i].main_circle ].x_center << "   " <<
            centers[i].circles[ centers[i].main_circle ].y_center << "   " <<
            centers[i].circles[ centers[i].main_circle ].radius << "   " <<
            centers[i].circles[ centers[i].main_circle ].n_votes << "   " <<
            centers[i].low_z << "   " << centers[i].up_z << endl;

        }

        result_file.close();

};

void saveCloud(vector<HoughCenters>* coordinates, double pos_z = 1.5, string file_path = "cloud.laz"){

                string format;
                unsigned pt = file_path.find_last_of(".");
                format = file_path.substr(pt+1);

                LASwriteOpener laswriteopener;
                laswriteopener.set_file_name(file_path.c_str());

                int format_macro;

                if(format == "las"){
                    format_macro = LAS_TOOLS_FORMAT_LAS;
                }else if(format == "laz"){
                    format_macro = LAS_TOOLS_FORMAT_LAZ;
                }else{
                    format_macro = LAS_TOOLS_FORMAT_TXT;
                }

                laswriteopener.set_format(format_macro);

                LASheader lasheader;
                lasheader.point_data_format = 0;
                lasheader.point_data_record_length = 20;
                LASpoint laspoint;
                laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, &lasheader);
                LASwriter* laswriter = laswriteopener.open(&lasheader);

                laspoint.set_z(pos_z);

                vector<HoughCenters>::iterator point;
                point = coordinates->begin();
                while (point != coordinates->end()){

                    vector<HoughCircle>::iterator c_point;
                    c_point = point->circles.begin();

                    while(c_point != point->circles.end()){

                        laspoint.set_x( (*c_point).x_center );
                        laspoint.set_y( (*c_point).y_center );
                        laswriter->write_point(&laspoint);
                        laswriter->update_inventory(&laspoint);

                        c_point++;
                    }
                    point++;
                }

                laswriter->update_header(&lasheader, TRUE);
                laswriter->close(TRUE);

                delete laswriter;

}

void printHelp(){

    cout <<
        "\n# /*** TLStools - las2rings ***/\n# /*** Command line arguments ***/\n\n"
        "# -i --input         : input file path\n"
        "# -o --output        : output file path (.txt)\n"
        "# -t --tree          : is single tree\n"
        "# -s --one-slice     : take only one slice\n"
        "# -l --lower         : slice's lower height\n"
        "# -u --upper         : slice's upper height\n"
        "# -p --pixel         : pixel size, in meters\n"
        "# -r --radius        : maximum radius to test\n"
        "# -d --density       : minimum density to consider on the Hough transform\n"
        "# -v --votes         : minimum votes count at the output\n"
        "# -O --output-cloud  : save a las/laz/txt cloud output\n"
        "# -? -h --help       : print help\n";

        exit(1);

}

/*************************************************************************************************************/
/** command line arguments **/

struct CommandLine {

    string file_path;            /* -i option, input file */
    string output_path;          /* -o option, output file */
    bool single_tree;            /* -t option, is single tree */
    bool one_slice;              /* -s option, take only one slice */
    string lower_slice;          /* -l option, lower height */
    string upper_slice;          /* -u option, upper height */
    double  pixel_size;          /* -p option, pixel size in meters */
    double  max_radius;          /* -r option, maximum radius to test */
    double  min_density;         /* -d option, minimum density to consider on the Hough transform */
    int  min_votes;              /* -v option, minimum votes count at the output */
    string output_las;           /* -O option, save a las/laz/txt output */
    bool help;                   /* -h/? option, get help */

} globalArgs;

static const char *optString = "i:o:tsl:u:p:r:d:v:O:h?";

static const struct option longOpts[] = {
    { "input", required_argument, NULL, 'i' },
    { "output", required_argument, NULL, 'o' },
    { "tree", no_argument, NULL, 't' },
    { "one-slice", no_argument, NULL, 's' },
    { "lower", required_argument, NULL, 'l' },
    { "upper", required_argument, NULL, 'u' },
    { "pixel", required_argument, NULL, 'p' },
    { "radius", required_argument, NULL, 'r' },
    { "density", required_argument, NULL, 'd' },
    { "votes", required_argument, NULL, 'v' },
    { "output-cloud", required_argument, NULL, 'O' },
    { "help", no_argument, NULL, 'h' },
    {NULL, no_argument, NULL, 0}
};

/*************************************************************************************************************/

int main(int argc, char *argv[])
{
    /** define global variables **/
    globalArgs.file_path = " ";
    globalArgs.help = false;
    globalArgs.lower_slice = "1.0";
    globalArgs.max_radius = 0.5;
    globalArgs.min_density = 0.1;
    globalArgs.min_votes = 3;
    globalArgs.one_slice = true;
    globalArgs.output_las = " ";
    globalArgs.output_path = " ";
    globalArgs.pixel_size = 0.025;
    globalArgs.single_tree = false;
    globalArgs.upper_slice = "2.0";

    int opt = 0;
	int longIndex = 0;

    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
            case 'i':
                globalArgs.file_path = std::string(optarg);
                break;

            case 'o':
                globalArgs.output_path = std::string(optarg);
                break;

            case 't':
                globalArgs.single_tree = true;
                break;

            case 's':
                globalArgs.one_slice = true;
                break;

            case 'l':
                globalArgs.lower_slice = std::string(optarg);
                break;

            case 'u':
                globalArgs.upper_slice = std::string(optarg);
                break;

            case 'p':
                globalArgs.pixel_size = atof(optarg);
                break;

            case 'r':
                globalArgs.max_radius = atof(optarg);
                break;

            case 'd':
                globalArgs.min_density = atof(optarg);
                break;

            case 'v':
                globalArgs.min_votes = atoi(optarg);
                break;

            case 'O':
                globalArgs.output_las = std::string(optarg);
                break;

            case 'h':
            case '?':
                globalArgs.help = true;
                break;

            default:
                break;
        }

        opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    }

    /** print help **/
    if(globalArgs.help){
        printHelp();
        return 0;
    }


    //globalArgs.file_path = "lcer.las";


    if(globalArgs.file_path == " "){
        cout << "\n# input file missing.\n";
        printHelp();
        return 0;
    }

    /** rename output **/
    if(globalArgs.output_path == " "){
        globalArgs.output_path = outputNameAppend(globalArgs.file_path);
    }

    if(globalArgs.output_las == " "){
        globalArgs.output_las = outputNameAppend(globalArgs.file_path, "_cloud.las");
    }

    cout << "# reading point cloud" << endl;
    Slice slc = getSlice(globalArgs.file_path, globalArgs.lower_slice, globalArgs.upper_slice);

    cout << "# rasterizing cloud's slice" << endl;
    Raster ras = getCounts(&slc, globalArgs.pixel_size);

    cout << "# extracting center candidates" << endl;
    vector<HoughCenters> hough = getCenters(&ras, globalArgs.max_radius, globalArgs.min_density, globalArgs.min_votes);

    cout << "# extracting precise center estimates" << endl;
    getPreciseCenters(hough);

    if(globalArgs.single_tree){
        HoughCircle vt;
        vt.n_votes = 0;
        for(unsigned i = 0; i < hough.size(); ++i){
            if(hough[i].circles[ hough[i].main_circle ].n_votes > vt.n_votes) vt = hough[i].circles[ hough[i].main_circle ];
        }
        cout << "# n of votes: " << vt.n_votes << ", radius: " << vt.radius << endl;
    }

    cout << "# writing cloud of center candidates: " << globalArgs.output_las << endl;
    saveCloud(&hough, (ras.max_z + ras.min_z)/2, globalArgs.output_las);

    cout << "# writing results: " << globalArgs.output_path << endl;
    saveReport(hough, globalArgs.output_path);

    cout << "# done" << endl;

    return 0;
}
