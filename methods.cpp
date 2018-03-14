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

//---------------------
//general methods
//---------------------


//output suffix
string outputNameAppend(string path, string suffix){

    unsigned last_dot = path.find_last_of(".");
    string new_path = path.substr(0, last_dot) + suffix;

    return new_path;
}

//check file
void checkInput(string filePath){

  LASreadOpener lasreadopener;
  lasreadopener.set_file_name(filePath.c_str());
  LASreader* lasreader = lasreadopener.open();
  lasreader->close();
  delete lasreader;

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

template <typename NumberType>
string number2string(NumberType number){

    stringstream ss;
    ss << number;

    string str;
    str = ss.str();

    return str;
}

float text2float(string number){

    stringstream ss;
    ss << number;

    float num;
    ss >> num;

    return num;
}


//---------------------
//plot and tree methods
//---------------------

//returns all points in between two heights and their x,y range
Slice getSlice(string file, string lower, string upper, float zfloor, bool clipTree, float xCenter, float yCenter, float centerRadius){

    float lpf = atof(lower.c_str()) + zfloor;
    float upf = atof(upper.c_str()) + zfloor;

    std::ostringstream ssl;
    ssl << lpf;
    lower = ssl.str();

    std::ostringstream ssu;
    ssu << upf;
    upper = ssu.str();

    //cout << "Z Range: " << lower << " : " << upper << endl;

    unsigned last_dot = file.find_last_of(".");
    string file_format = file.substr(last_dot+1, file.length());

    Slice slc;
    slc.slice = {};

    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(file.c_str());

    char* MY_other_argv[4];

    MY_other_argv[0] = (char*)malloc(1000); strcpy(MY_other_argv[0], " ");
    MY_other_argv[1] = (char*)malloc(1000); strcpy(MY_other_argv[1], "-keep_z");
    MY_other_argv[2] = (char*)malloc(1000); strcpy(MY_other_argv[2], lower.c_str() );
    MY_other_argv[3] = (char*)malloc(1000); strcpy(MY_other_argv[3], upper.c_str() );

    lasreadopener.parse(4, MY_other_argv);
    LASreader* lasreader = lasreadopener.open();

    unsigned i = 0;
    double min_x, max_x, min_y, max_y, min_z, max_z;

    double x;
    double y;
    double z;
    double xyDist;

    while(lasreader->read_point()){

        x = lasreader->get_x();
        y = lasreader->get_y();
        z = lasreader->get_z();

        if(clipTree){
            xyDist = sqrt( pow(x - xCenter, 2) + pow(y - yCenter, 2) );
            if(xyDist > centerRadius) continue;
        }

        if(i == 0){
           min_x = max_x = x;
           min_y = max_y = y;
           min_z = max_z = z;
        }

        slc.slice.push_back( {x,y,z} );

        if(x > max_x) max_x = x;
        if(y > max_y) max_y = y;
        if(x < min_x) min_x = x;
        if(y < min_y) min_y = y;
        if(z < min_z) min_z = z;
        if(z > max_z) max_z = z;

        i++;

    }

    lasreader->close();
    delete lasreader;

    slc.dims = {min_x, max_x, min_y, max_y, min_z, max_z};

    return slc;

}

//returns a matrix of counts and maximum count value
Raster getCounts(Slice* slice , float pixel_size, double x_mid, double y_mid, double d_mid){

    Raster ras;
    ras.min_x = (d_mid > 0) ? (x_mid - d_mid) : (*slice).dims[0];
    ras.max_x = (d_mid > 0) ? (x_mid + d_mid) : (*slice).dims[1];
    ras.min_y = (d_mid > 0) ? (y_mid - d_mid) : (*slice).dims[2];
    ras.max_y = (d_mid > 0) ? (y_mid + d_mid) : (*slice).dims[3];
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
    double x,y;
    int xCell;
    int yCell;

    for(unsigned int i = 0; i < (*slice).slice.size(); i++){

        x = (*slice).slice[i][0];
        y = (*slice).slice[i][1];

        double dist = sqrt( pow(x - x_mid,2) + pow(y - y_mid,2) );
        if(d_mid > 0 && dist > d_mid) continue;

        xCell = abs( floor( (x - ras.min_x) / pixel_size ) );
        yCell = abs( floor( (y - ras.min_y) / pixel_size ) );

        if(xCell >= xn) xCell = xn-1;
        if(yCell >= yn) yCell = yn-1;

        matrix[xCell][yCell] += 1;

        if(matrix[xCell][yCell] > max_votes) max_votes = matrix[xCell][yCell];
    }

    ras.matrix = matrix;
    ras.max_count = max_votes;

    return ras;

}

//create a rasterized circle
vint rasterCircle(float radius, float pixel_size, double cx, double cy, double mx, double my){

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

vector<HoughCenters> getCenters(Raster* raster, float max_radius, float min_den, int min_votes){

    //count raster properties (&raster)
    unsigned int x_len = raster->x_dim;
    unsigned int y_len = raster->y_dim;
    int min_count = ceil( (raster->max_count)*min_den );

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

void saveReport(vector<HoughCenters>& centers, string file_path){

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

void saveCloud(vector<HoughCenters>* coordinates, string file_path){

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

                //laspoint.set_z(pos_z);

                vector<HoughCenters>::iterator point;
                point = coordinates->begin();
                while (point != coordinates->end()){

                    vector<HoughCircle>::iterator c_point;
                    c_point = point->circles.begin();
                    laspoint.set_z( (point->low_z + point->up_z)/2 );

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

void saveStemsOnly(vector<vector<StemSegment>>& stemsList, float pixel, string input_path, string file_path, string cloud_path){

    //TXT writer
    ofstream result_file(file_path);

    result_file << "tree   x   y   rad   votes   z_min   z_max   points" << endl;


    //write files
    int counter = 0;

    for(auto i = stemsList.begin(); i != stemsList.end(); ++i){

        counter++;

        for(auto j = i->begin(); j != i->end(); ++j){
           result_file <<
                counter << "   " <<
                j->model_circle.x_center << "   " <<
                j->model_circle.y_center << "   " <<
                j->model_circle.radius << "   " <<
                j->model_circle.n_votes << "   " <<
                j->z_min << "   " <<
                j->z_max << "   " <<
                j->n_points << "   "
            << endl;

        }
    }

    result_file.close();


    //LAS read/write
    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(input_path.c_str());
    LASreader* lasreader = lasreadopener.open();

    LASwriteOpener laswriteopener;
    laswriteopener.set_file_name(cloud_path.c_str());
    laswriteopener.set_format(LAS_TOOLS_FORMAT_LAZ);

    LASheader lasheader;
    lasheader.point_data_format = 0;
    lasheader.point_data_record_length = 20;
    LASpoint laspoint;
    laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, &lasheader);
    LASwriter* laswriter = laswriteopener.open(&lasheader);

    while (lasreader->read_point()){

        bool breakList = false;

        float tempX = lasreader->get_x();
        float tempY = lasreader->get_y();
        float tempZ = lasreader->get_z();

        for(auto i = stemsList.begin(); i != stemsList.end(); ++i){
            for(auto j = i->begin(); j != i->end(); ++j){

                if(tempZ >= j->z_min && tempZ <= j->z_max){

                    float pointDist = sqrt( pow(tempX - j->model_circle.x_center, 2) + pow(tempY - j->model_circle.y_center, 2) );

                    if(pointDist <= j->model_circle.radius + pixel){
                       laspoint.set_x( tempX );
                       laspoint.set_y( tempY );
                       laspoint.set_z( tempZ );

                       //cout << tempX << ", " << tempY << ", " << tempZ << endl;

                       laswriter->write_point(&laspoint);
                       laswriter->update_inventory(&laspoint);

                       breakList = true;
                       break;
                    }

                }

            }

            if(breakList){
                break;
            }

        }

    }

    laswriter->close();
    delete laswriter;

    lasreader->close();
    delete lasreader;

}

vector<vector<HoughCircle*>> isolateSingleTrees(vector<HoughCenters>& roughTreeMap, float maxDist, int minLayers){

    vector<vector<HoughCircle*>> circleClusters;
    HoughCircle* temp = &roughTreeMap[0].circles[ roughTreeMap[0].main_circle ];

    circleClusters.push_back( { temp } );

    for(int i = 1; i < roughTreeMap.size(); ++i){

        temp = &roughTreeMap[i].circles[ roughTreeMap[i].main_circle ];

        bool broken = false;
        for(vector<vector<HoughCircle*>>::iterator j = circleClusters.begin(); j != circleClusters.end(); j++ ){
            for(vector<HoughCircle*>::iterator k = j->begin(); k != j->end(); k++){

                float euclidean = sqrt( pow(temp->x_center - (**k).x_center , 2) + pow(temp->y_center - (**k).y_center , 2) );
                if(euclidean <= maxDist){
                    j->push_back( temp );
                    broken = true;
                    break;
                }

            }

            if(broken){
                break;
            }
        }

        if(!broken){
            circleClusters.push_back( { temp } );
        }

    }

    for(int i = circleClusters.size() - 1; i >= 0; --i){
        if( circleClusters[i].size() < minLayers){
            circleClusters.erase( circleClusters.begin() + i );
        }
    }

    return circleClusters;
}

////////////////////////
// single tree functions
////////////////////////

vector<Slice> sliceList(string file, CloudStats& props, float z_interval, bool clipTree, float xCenter, float yCenter, float centerRadius){

    vector<double> heights = {};
    Slice s1;
    s1.slice = {};
    int n_slices = ceil( (props.z_max - props.z_min) / z_interval );
    vector<Slice> slices(n_slices, s1);

    int si = 0;
    for(double i = props.z_min; i < props.z_max + z_interval; i+=z_interval){
        heights.push_back(i);
        if(i < props.z_max && !clipTree){

          slices[si++].dims = {props.x_max, props.x_min, props.y_max, props.y_min, i, i+z_interval};

        }else if(i < props.z_max && clipTree){

          slices[si++].dims = {xCenter, xCenter, yCenter, yCenter, i, i+z_interval};

        }
    }

    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(file.c_str());
    LASreader* lasreader = lasreadopener.open();

    double x, y, z;
    double xyDist;


    while(lasreader->read_point()){

        x = lasreader->get_x();
        y = lasreader->get_y();
        z = lasreader->get_z();

        if(clipTree){
            xyDist = sqrt( pow(x - xCenter, 2) + pow(y - yCenter, 2) );
            if(xyDist > centerRadius) continue;
        }

        si = floor( (z - props.z_min) / z_interval );
        if(si < 0) si = 0;
        if(si > (n_slices - 1)) si = (n_slices - 1);
        slices[si].slice.push_back( {x,y,z} );

        if(x < slices[si].dims[0]) slices[si].dims[0] = x;
        if(x > slices[si].dims[1]) slices[si].dims[1] = x;
        if(y < slices[si].dims[2]) slices[si].dims[2] = y;
        if(y > slices[si].dims[3]) slices[si].dims[3] = y;

    }

    lasreader->close();
    delete lasreader;

    vector<Slice> returnSlices;
    for(std::vector<Slice>::iterator it = slices.begin(); it != slices.end(); ++it){
        if(it->slice.size() < 5) continue;
        returnSlices.push_back( (*it) );
    }

    return returnSlices;

}

int getMainEstimate(vector<HoughCenters>& circlesList){

    int nv = 0;
    int found;
    for( unsigned i=0; i < circlesList.size(); ++i ){
        if(circlesList[i].getMainCircle().n_votes > nv){
            nv = circlesList[i].getMainCircle().n_votes;
            found = i;
        }
    };

    return found;

}

StemSegment baselineStats(CloudStats& stats, CommandLine global, bool clipTree, float xCenter, float yCenter, float centerRadius){
/*
    double h1 = atof(global.lower_slice.c_str());
    double h2 = atof(global.upper_slice.c_str());

    //convert double to string
    std::ostringstream ssh1, ssh2;
    ssh1 << h1;
    ssh2 << h2;
    std::string sh1(ssh1.str());
    std::string sh2(ssh2.str());
*/
    Slice base =  getSlice(global.file_path, global.lower_slice, global.upper_slice, stats.z_min, clipTree, xCenter, yCenter, centerRadius);

    Raster raster = getCounts(&base, global.pixel_size);

    vector<HoughCenters> circles = getCenters(&raster, global.max_radius, global.min_density, global.min_votes);
    getPreciseCenters(circles);

    StemSegment stem_base;

    int mc = getMainEstimate(circles);
    stem_base.model_circle = circles[mc].getMainCircle();
    stem_base.n_points = base.slice.size();
    stem_base.z_max = raster.max_z;
    stem_base.z_min = raster.min_z;

    return stem_base;

};

StemSegment getSegment(Slice& slc, CommandLine global, float xm, float ym, float dm){

    StemSegment stSeg;

    if( slc.slice.size() < 5){
        return stSeg;
    }

    Raster segRas = getCounts(&slc, global.pixel_size, xm, ym, dm);

    vector<HoughCenters> segHough = getCenters(&segRas, dm, global.min_density, global.min_votes);

    getPreciseCenters(segHough);

    if(segHough.size() == 0){
        return stSeg;
    }

    int nv = 0;
    int found;
    for( unsigned i=0; i < segHough.size(); ++i ){
        if(segHough[i].getMainCircle().n_votes > nv){
            nv = segHough[i].getMainCircle().n_votes;
            found = i;
        }
    };

    stSeg.model_circle = segHough[found].getMainCircle();
    stSeg.n_points = slc.slice.size();
    stSeg.z_min = segRas.min_z;
    stSeg.z_max = segRas.max_z;

    return stSeg;

};

void saveStemReport(vector<StemSegment>& sections, string file_path){

        ofstream stem_file(file_path);

        stem_file << "x_center   y_center   radius   votes   z_min   z_max" << endl;

        for(unsigned int i = 0; i < sections.size() ; ++i){

            stem_file <<
                sections[i].model_circle.x_center << "   " <<
                sections[i].model_circle.y_center << "   " <<
                sections[i].model_circle.radius << "   " <<
                sections[i].model_circle.n_votes << "   " <<
                sections[i].z_min << "   " <<
                sections[i].z_max
            << endl;

        }

        stem_file.close();

}

void saveStemCloud(vector<StemSegment>& stem, vector<Slice>& tree, double pixel_size,string file_path){

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

                for(unsigned i = 0; i < stem.size(); ++i){
                    double dist = stem[i].model_circle.radius + pixel_size;

                    for(unsigned j = 0; j < tree[i].slice.size(); ++j){

                        double x = tree[i].slice[j][0];
                        double y = tree[i].slice[j][1];
                        double z = tree[i].slice[j][2];

                        double d = sqrt( pow(x - stem[i].model_circle.x_center ,2) + pow(y - stem[i].model_circle.y_center ,2) );

                        if(d <= dist){
                            laspoint.set_x( x );
                            laspoint.set_y( y );
                            laspoint.set_z( z );

                            laswriter->write_point(&laspoint);
                            laswriter->update_inventory(&laspoint);
                        }
                    }
                }


                laswriter->update_header(&lasheader, TRUE);
                laswriter->close(TRUE);

                delete laswriter;

}

vector<StemSegment> stemPoints(StemSegment& base, vector<Slice>& pieces, CommandLine global){

    float xt = base.model_circle.x_center;
    float yt = base.model_circle.y_center;
    float dt = base.model_circle.radius + global.pixel_size*4;
    vector<StemSegment> stem_sections = {};

    for(unsigned i = 0; i < pieces.size(); ++i){

        //cout << "stem slice " << i << endl;

        Slice* temp = &pieces[i];

        Raster ras = getCounts(temp, global.pixel_size);

        vector<HoughCenters> hough = getCenters(&ras, dt, global.min_density, global.min_votes);

        getPreciseCenters(hough);

        if(hough.size() < 1) continue;

        int nCircle = 0;
        int nVotes  = -1;
        for(int j = 0; j < hough.size(); ++j){
            HoughCenters* it = &hough[j];
            int nv = it->circles[it->main_circle].n_votes;
            if(nv > nVotes){
                nVotes  = nv;
                nCircle = j;
            }
        }

        int id = hough[nCircle].main_circle;

        if(ras.max_z > base.z_max){
            dt = hough[nCircle].circles[id].radius + global.pixel_size*1;
        }

        StemSegment mainBolePiece;
        mainBolePiece.model_circle = hough[nCircle].circles[id];
        mainBolePiece.n_points     = temp->slice.size();
        mainBolePiece.z_max        = ras.max_z;
        mainBolePiece.z_min        = ras.min_z;

        stem_sections.push_back(mainBolePiece);
    }

    return stem_sections;

}
