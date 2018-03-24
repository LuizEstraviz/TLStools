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
#include "methods.hpp"

using namespace std;

/*************************************************************************************************************/
CommandLine globalArgs;

static const char *optString = "i:o:l:u:p:r:d:v:O:z:h?";

static const struct option longOpts[] = {
    { "input", required_argument, NULL, 'i' },
    //{ "output", required_argument, NULL, 0 },
    { "lower", required_argument, NULL, 'l' },
    { "upper", required_argument, NULL, 'u' },
    { "pixel", required_argument, NULL, 'p' },
    { "radius", required_argument, NULL, 'r' },
    { "density", required_argument, NULL, 'd' },
    { "votes", required_argument, NULL, 'v' },
    //{ "output-cloud", required_argument, NULL, 0 },
    { "help", no_argument, NULL, 'h' },
    { "zheight", required_argument, NULL, 'z' },
    //{ "stack-stats", required_argument, NULL, 0},
    //{ "stack-cloud", required_argument, NULL, 0},
    { "olazdir", required_argument, NULL, 'O'},
    { "otxtdir", required_argument, NULL, 'o'},
    {NULL, no_argument, NULL, 0}
};

/*************************************************************************************************************/

/*** general functions ***/

void printHelp(){

    cout <<
        "\n# /*** TLStrees ***/\n# /*** Command line arguments ***/\n\n"
        "# -i or --input         : input file path\n"
        "# -o or --otxtdir       : output directory for the report files (defaults to current directory)\n"
        "# -O or --olazdir       : output directory for the point cloud files (defaults to current directory)\n"
        "# -l or --lower         : slice's lower height (default = 1.0 m)\n"
        "# -u or --upper         : slice's upper height (default = 3.0 m)\n"
        "# -z or --zheight       : height interval to search for stem segments (default = 0.5 m)\n"
        "# -p or --pixel         : pixel size, in meters (default = 0.025 m)\n"
        "# -r or --radius        : maximum radius to test (default = 0.25 m)\n"
        "# -d or --density       : minimum density to consider on the Hough transform, from 0 to 1 (default = 0.1)\n"
        "# -v or --votes         : minimum votes count at the output (default = 3 votes per pixel)\n"
        "# -? or -h or --help    : print this help\n"
        //"# --stack-stats         : output file path for the tree center statistics (.txt)\n"
        //"# --stack-cloud         : output file path for the trees center layer stack (.laz/.las/.txt)\n"
    << endl;

        exit(1);

}

/*** batch processing ***/
//statistics per tree
vector<PlotCloudTreeSlice> getTlsMetrics(vector<HoughCenters>& centers){

    vector<PlotCloudTreeSlice> tlsMetrics;
    PlotCloudTreeSlice metric;

    for(unsigned int i = 0; i < centers.size() ; ++i){

        metric.n_tree = i+1;
        metric.n_slice = 0;
        metric.x_average = centers[i].avg_x;
        metric.y_average = centers[i].avg_y;
        metric.x_main = centers[i].circles[ centers[i].main_circle ].x_center;
        metric.y_main = centers[i].circles[ centers[i].main_circle ].y_center;
        metric.radius = centers[i].circles[ centers[i].main_circle ].radius;
        metric.votes = centers[i].circles[ centers[i].main_circle ].n_votes;
        metric.z_min = centers[i].low_z;
        metric.z_max = centers[i].up_z;

        tlsMetrics.push_back(metric);
    }

    return tlsMetrics;

}
/**/
//plot-wise
void plotProcess(CommandLine global){

    vector<HoughCenters> treeMap;
    Raster ras;

    cout << "# getting cloud statistics" << endl;
    CloudStats cstats = getStats(global.file_path);

    float hLow = text2float(global.lower_slice);
    float hUp = text2float(global.upper_slice);

    cout << "# mapping tree positions" << endl;
    for(float h = hLow; h <= (hUp - global.height_interval); h += global.height_interval){

        cout << "## layer " << 1 + ( h - hLow ) / global.height_interval << endl;

        std::stringstream ss, ss2;
        ss << h;
        ss2 << (h+global.height_interval);
        string ht;
        string ht2;
        ss >> ht;
        ss2 >> ht2;

        //cout << "# reading point cloud" << endl;
        Slice slc = getSlice(global.file_path, ht, ht2, cstats.z_min);

        //cout << "# rasterizing cloud's slice" << endl;
        ras = getCounts(&slc, global.pixel_size);

        //cout << "# extracting center candidates" << endl;
        vector<HoughCenters> hough = getCenters(&ras, global.max_radius, global.min_density, global.min_votes);

        //cout << "# extracting center estimates" << endl;
        getPreciseCenters(hough);

        //treeMap.resize(treeMap.size() + hough.size());
        treeMap.insert(treeMap.end(), hough.begin(), hough.end());

    }

    cout << "# writing cloud of center candidates: " << global.output_stack << endl;
    saveCloud(&treeMap, global.output_stack);

    cout << "# writing layer stack results: " << global.output_stack_coordinates << endl;
    saveReport(treeMap, global.output_stack_coordinates);

    cout << "# extracting main coordinate per tree" << endl;
    int nLayers = 0.75 * (hUp - hLow) / global.height_interval;
    vector<vector<HoughCircle*>> singleTreeMap = isolateSingleTrees(treeMap, global.max_radius * 2, nLayers);

    cout << "# getting trunk statistics" << endl;
    vector<vector<StemSegment>> trees;
    for(int i = 0; i < singleTreeMap.size(); ++i){
        cout << "\ntree " << i+1 << " of " << singleTreeMap.size() << endl;

        vector<HoughCircle*> temp = singleTreeMap[i];

        float xBase = temp[0]->x_center;
        float yBase = temp[0]->y_center;
        cout << "center: " << xBase << " , " << yBase << endl;

        StemSegment base = baselineStats(cstats, global, true, xBase, yBase, 1.2);

        vector<Slice> pieces = sliceList(global.file_path, cstats, global.height_interval, true, xBase, yBase, global.max_radius*3);
        cout << "segments: " << pieces.size() << endl;

        vector<StemSegment> bole = stemPoints(base, pieces, global);
        trees.push_back(bole);
    }

    cout << endl;

    cout << "# writing tree-wise statistics to " << global.output_path << " and stem point cloud to " << global.output_las << endl;
    saveStemsOnly(trees, global.pixel_size, global.file_path, global.output_path, global.output_las, global.max_radius);

    cout << "# done" << endl;

}

int main(int argc, char *argv[])
{

    /** define global variables **/
    globalArgs.file_path = " ";
    globalArgs.help = false;
    globalArgs.lower_slice = "1.0";
    globalArgs.max_radius = 0.25;
    globalArgs.min_density = 0.1;
    globalArgs.min_votes = 3;
    globalArgs.output_las = " ";
    globalArgs.output_path = " ";
    globalArgs.pixel_size = 0.025;
    globalArgs.upper_slice = "3.0";
    globalArgs.height_interval = 0.5;
    globalArgs.output_stack = " ";
    globalArgs.output_stack_coordinates = " ";
    globalArgs.clouds_directory = "";
    globalArgs.reports_directory = "";

    int opt = 0;
	int longIndex = 0;

    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
            case 'i':
                globalArgs.file_path = std::string(optarg);
                break;

            case 'o':
                globalArgs.reports_directory = std::string(optarg);
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
                globalArgs.clouds_directory = std::string(optarg);
                break;

            case 'z':
                globalArgs.height_interval = atof(optarg);
                break;

            case 'h':
            case '?':
                globalArgs.help = true;
                break;

            /*case 0:

                if( strcmp( "stack-stats", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.output_stack_coordinates = std::string(optarg);

                }else if(strcmp( "stack-cloud", longOpts[longIndex].name ) == 0){
                    globalArgs.output_stack = std::string(optarg);

                }else if(strcmp( "output", longOpts[longIndex].name ) == 0){
                    globalArgs.output_path = std::string(optarg);

                }else if(strcmp( "output-cloud", longOpts[longIndex].name ) == 0){
                    globalArgs.output_las = std::string(optarg);

                }

                break;
            */

            default:
                break;
        }

        opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    }


    /*** parse command line options ***/

    // print helpvector<HoughCenters>&
    if(globalArgs.help){
        printHelp();
        return 0;
    }
/*
    globalArgs.file_path = "sample_data/square.las";
*/
    if(globalArgs.file_path == " "){
        cout << "\n# input file (-i) missing." << endl;
        //printHelp();
        return 0;
    }

    if( (text2float(globalArgs.upper_slice) - text2float(globalArgs.lower_slice)) < globalArgs.height_interval ){
        cout << "\n# the height interval (-z) cannot be larger than the distance between lower (-l) and upper (-u) limits!" << endl;
        //printHelp();
        return 0;
    }

    //parse names
    int lastSlash = globalArgs.file_path.find_last_of("/");
    string shortFileName = (lastSlash == -1) ? globalArgs.file_path : globalArgs.file_path.substr(lastSlash, globalArgs.file_path.length());

    string lastChar = globalArgs.reports_directory.substr( globalArgs.reports_directory.length() );
    if(lastChar != "/" && globalArgs.reports_directory != ""){
        globalArgs.reports_directory += "/";
    }

    lastChar = globalArgs.clouds_directory.substr( globalArgs.clouds_directory.length() );
    if(lastChar != "/" && globalArgs.clouds_directory != ""){
        globalArgs.clouds_directory += "/";
    }


    // rename output
    if(globalArgs.output_path == " "){
        globalArgs.output_path = outputNameAppend(shortFileName);
    }

    if(globalArgs.output_las == " "){
        globalArgs.output_las = outputNameAppend(shortFileName, "_trees.laz");
    }

    if(globalArgs.output_stack == " "){
        globalArgs.output_stack = outputNameAppend(shortFileName, "_segmt.laz");
    }

    if(globalArgs.output_stack_coordinates == " "){
        globalArgs.output_stack_coordinates = outputNameAppend(shortFileName, "_segmt.txt");
    }

    globalArgs.output_path              = globalArgs.reports_directory + globalArgs.output_path;
    globalArgs.output_las               = globalArgs.clouds_directory  + globalArgs.output_las;
    globalArgs.output_stack             = globalArgs.clouds_directory  + globalArgs.output_stack;
    globalArgs.output_stack_coordinates = globalArgs.reports_directory + globalArgs.output_stack_coordinates;

    /*** process ***/
    checkInput(globalArgs.file_path);
    plotProcess(globalArgs);

    return 0;
}
