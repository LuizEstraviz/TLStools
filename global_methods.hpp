#ifndef GLOBAL_METHODS_HPP_INCLUDED
#define GLOBAL_METHODS_HPP_INCLUDED

#include <string>
#include <vector>
#include "classes.hpp"
#include "global_methods.hpp"

using namespace std;

//output suffix
string outputNameAppend(string path, string suffix = "_result.txt");

//get cloud statistics
CloudStats getStats(string file);

//calculate absolute center coordinate based on pixel position
vector<double> absCenter(int x, int y, double min_x, double min_y, float step);

//calculate pixel based on absolute coordinate
vector<int> pixPosition(double x, double y, double min_x, double min_y, float step);

#endif // GLOBAL_METHODS_HPP_INCLUDED
