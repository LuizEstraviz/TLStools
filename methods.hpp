#ifndef METHODS_HPP_INCLUDED
#define METHODS_HPP_INCLUDED

#include <vector>
#include <string>
#include "classes.hpp"

using namespace std;

//output suffix
string outputNameAppend(string path, string suffix = "_result.txt");

//get cloud statistics
CloudStats getStats(string file);

//calculate absolute center coordinate based on pixel position
vector<double> absCenter(int x, int y, double min_x, double min_y, float step);

//calculate pixel based on absolute coordinate
vector<int> pixPosition(double x, double y, double min_x, double min_y, float step);

//returns all points in between two heights and their x,y range
Slice getSlice(string file, string lower = "1.0", string upper = "2.0", float zfloor = 0);

//returns a matrix of counts and maximum count value
Raster getCounts(Slice* slice , float pixel_size = 0.025, double x_mid = 0, double y_mid = 0, double d_mid = -1);

//create a rasterized circle
vint rasterCircle(float radius, float pixel_size = 0.025, double cx = 0, double cy = 0, double mx = 0, double my = 0);

vector<HoughCenters> getCenters(Raster* raster, float max_radius = 0.5, float min_den = 0.1, int min_votes = 2);

void getPreciseCenters(vector<HoughCenters>& circles);

void saveReport(vector<HoughCenters>& centers, string file_path = "result.txt");

void saveCloud(vector<HoughCenters>* coordinates, string file_path = "cloud.laz");


vector<Slice> sliceList(string file, CloudStats& props, float z_interval = 0.5);

int getMainEstimate(vector<HoughCenters>& circlesList);

StemSegment baselineStats(CloudStats& stats, CommandLine global);

StemSegment getSegment(Slice& slc, CommandLine global, float xm, float ym, float dm);

void saveStemReport(vector<StemSegment>& sections, string file_path = "stem.txt");

void saveStemCloud(vector<StemSegment>& stem, vector<Slice>& tree, double pixel_size=0.025 ,string file_path = "stem_cloud.laz");


#endif // METHODS_HPP_INCLUDED
