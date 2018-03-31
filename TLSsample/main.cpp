#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <sstream>
#include <getopt.h>
#include <random>
#include "lasreader.hpp"
#include "laswriter.hpp"

using namespace std;

/*************************************************************************************************************/
struct CommandLine {

    string file_path;
    string suffix;
    float voxel_size;
    bool help;
    float random_proportion;

} globalArgs;

static const char *optString = "i:o:v:r:h?";

static const struct option longOpts[] = {
    { "input", required_argument, NULL, 'i' },
    { "voxel", required_argument, NULL, 'v' },
    { "help", no_argument, NULL, 'h' },
    { "odix", required_argument, NULL, 'o'},
    { "random", required_argument, NULL, 'r'},
    {NULL, no_argument, NULL, 0}
};

/*************************************************************************************************************/

/*** general functions ***/

void printHelp(){

    cout <<
        "\n# /*** TLSsample ***/\n# /*** Command line arguments ***/\n\n"
        "# -i or --input         : input file path\n"
        "# -o or --odix          : output suffix for the sampled point cloud (default = sample)\n"
        "# -v or --voxel         : (default method) voxel side length (default = 0.05 m)\n"
        "# -r or --random        : randomly sample a proportion (0-1] of the input cloud.\n"
        "                         *this method is applied only if the argument is explicitly given. Otherwise, voxel sampling is performed.\n"
        "# -? or -h or --help    : print this help\n"
        //random
    << endl;

        exit(1);

}

/*************************************************************************************************************/


class minMax{
    public:
        float x[2];
        float y[2];
        float z[2];
};

template <typename NumberType>
string number2string(NumberType number){

    stringstream ss;
    ss << number;

    string str;
    str = ss.str();

    return str;
}

string renameSuffix(string path, string suffix = "trim"){

    unsigned last_dot = path.find_last_of(".");
    string new_path = path.substr(0, last_dot) + "_" + suffix + path.substr(last_dot, path.length());

    if(new_path[0] == '/') new_path = "." + new_path;

    return new_path;
}

minMax getLimits(string file, float voxelSize = 0.05){

  LASreadOpener lasreadopener;
  lasreadopener.set_file_name(file.c_str());
  LASreader* lasreader = lasreadopener.open();

  minMax cloudRange;

  unsigned last_dot = file.find_last_of(".");
  string file_format = file.substr(last_dot+1, file.length());

  if(file_format == "las" || file_format == "laz"){

    cloudRange.x[0] = lasreader->get_min_x() - voxelSize;
    cloudRange.x[1] = lasreader->get_max_x() + voxelSize;

    cloudRange.y[0] = lasreader->get_min_y() - voxelSize;
    cloudRange.y[1] = lasreader->get_max_y() + voxelSize;

    cloudRange.z[0] = lasreader->get_min_z() - voxelSize;
    cloudRange.z[1] = lasreader->get_max_z() + voxelSize;

  } else {

    unsigned i = 1;
    float x, y, z, min_x, max_x, min_y, max_y, min_z, max_z;

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

    cloudRange.x[0] = min_x - voxelSize;
    cloudRange.x[1] = max_x + voxelSize;

    cloudRange.y[0] = min_y - voxelSize;
    cloudRange.y[1] = max_y + voxelSize;

    cloudRange.z[0] = min_z - voxelSize;
    cloudRange.z[1] = max_z + voxelSize;

    }

    lasreader->close();
    delete lasreader;

    return cloudRange;

};

string getVoxelString(float x, float y, float z, minMax& range, float voxelSize = 0.05){

  int fixedFactor = 100000;

  long long int nx = fixedFactor * ( floor( (x - range.x[0]) / voxelSize ) / voxelSize + range.x[0] );
  long long int ny = fixedFactor * ( floor( (y - range.y[0]) / voxelSize ) / voxelSize + range.y[0] );
  long long int nz = fixedFactor * ( floor( (z - range.z[0]) / voxelSize ) / voxelSize + range.z[0] );

  string key = number2string(nx) + "_" + number2string(ny) + "_" + number2string(nz);

  return key;

}

void voxelSample(string path, string suffix, float voxelSize = 0.05){

    minMax cloudRange = getLimits(path, voxelSize);
    string sampled_file = renameSuffix(path, suffix);
    unordered_map<string, short int> ledger;

    #ifdef DEBUGMODE
        cout << "output: " << sampled_file << endl;
        cout << "x: " << cloudRange.x[0] << " -> " << cloudRange.x[1] << endl;
        cout << "y: " << cloudRange.y[0] << " -> " << cloudRange.y[1] << endl;
        cout << "z: " << cloudRange.z[0] << " -> " << cloudRange.z[1] << endl;
    #endif // DEBUGMODE

    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(path.c_str());
    LASreader* lasreader = lasreadopener.open();

    LASwriteOpener laswriteopener;
    laswriteopener.set_file_name( sampled_file.c_str() );
    LASwriter* laswriter = laswriteopener.open(&lasreader->header);

    cout << "## writing " << sampled_file << endl;

    while(lasreader->read_point()){

        float x = lasreader->get_x();
        float y = lasreader->get_y();
        float z = lasreader->get_z();

        string key = getVoxelString(x, y, z, cloudRange, voxelSize);

        if( ledger.find(key) == ledger.end() ){
            ledger[key] = 1;
            laswriter->write_point(&lasreader->point);
            laswriter->update_inventory(&lasreader->point);
        }

    }

    laswriter->update_header(&lasreader->header, TRUE);
    laswriter->close();
    delete laswriter;

    lasreader->close();
    delete lasreader;
}

void randomSample(string path, string suffix, default_random_engine seed, float proportion = 0.1){

  uniform_real_distribution<float> distribution(0,1);

  string sampled_file = renameSuffix(path, suffix);

    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(path.c_str());
    LASreader* lasreader = lasreadopener.open();

    LASwriteOpener laswriteopener;
    laswriteopener.set_file_name( sampled_file.c_str() );
    LASwriter* laswriter = laswriteopener.open(&lasreader->header);

    cout << "## writing " << sampled_file << endl;

    while(lasreader->read_point()){

        float draw = distribution(seed);

        if( draw <= proportion ){
            laswriter->write_point(&lasreader->point);
            laswriter->update_inventory(&lasreader->point);
        }

    }

    laswriter->update_header(&lasreader->header, TRUE);
    laswriter->close();
    delete laswriter;

    lasreader->close();
    delete lasreader;

}

int main(int argc, char *argv[])
{
    globalArgs.help = false;
    globalArgs.suffix = "sample";
    globalArgs.voxel_size = 0.05;
    globalArgs.file_path = " ";
    globalArgs.random_proportion = 0;

    int opt = 0;
	int longIndex = 0;

    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
            case 'i':
                globalArgs.file_path = std::string(optarg);
                break;

            case 'o':
                globalArgs.suffix= std::string(optarg);
                break;

            case 'v':
                globalArgs.voxel_size = atof(optarg);
                break;

            case 'r':
                globalArgs.random_proportion = atof(optarg);
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

    /*** parse command line options ***/

    //globalArgs.file_path = "../sample_data/square.las";

    if(globalArgs.help){
        printHelp();
        return 0;
    }

    if(globalArgs.file_path == " "){
        cout << "\n# input file (-i) missing." << endl;
        return 0;
    }

    if(globalArgs.voxel_size <= 0){
        cout << "\n# the voxel size (-v) must be larger than 0." << endl;
        return 0;
    }

    if(globalArgs.random_proportion < 0 || globalArgs.random_proportion > 1){
        cout << "\n# the sampling proportion (-r) must be in the (0-1] range." << endl;
        return 0;
    }

    if(globalArgs.random_proportion == 0){
        cout << "## applying voxel sampling" << endl;

        voxelSample(globalArgs.file_path, globalArgs.suffix, globalArgs.voxel_size);

    }else{
        cout << "## applying random sampling" << endl;

        default_random_engine generator;
        randomSample(globalArgs.file_path, globalArgs.suffix, generator, globalArgs.random_proportion);
    }

    cout << "## done" << endl;

    return 0;
}
