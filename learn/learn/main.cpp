#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <iterator>

using namespace std;

string coor;
vector<double> x, y, z;
double former_string;
int i = 0;
int tx, ty, tz;

/*
bool  inPixel(double val, double x1, double x2, double y1, double y2){
    return  (val > 0 ? value  : -value) <= numeric_limits<RangeType>::max();
}
*/


int main()
{
  cout.precision(12);

  ifstream tab("lcer.txt");

  while( tab >> coor ) //getline (tab, ln) )
  {
    ++i ;

    if( i == 1 || i == tx){
        tx = i+3 ;
        former_string = strtod(coor.c_str(), NULL) ;
        x.push_back( former_string ) ;
    }

    else if(i == 2 || i == ty){
        ty = i+3 ;
        former_string = strtod(coor.c_str(), NULL) ;
        y.push_back( former_string ) ;

    }

    else if(i == 3 || i == tz){
        tz = i+3 ;
        former_string = strtod(coor.c_str(), NULL) ;
        z.push_back( former_string ) ;

    }

    else continue ;

  }

   tab.close() ;


   double H_low = 1, H_up = 2;
   vector<double> x_slice, y_slice, z_slice;

   for(unsigned i = 0 ; i != x.size() ; i++){
    if(z[i] >= H_low && z[i] <= H_up ){
        x_slice.push_back(x[i]);
        y_slice.push_back(y[i]);
        z_slice.push_back(z[i]);
    }
   }

   double min_x = *min_element(begin(x_slice), end(x_slice));
   double max_x = *max_element(begin(x_slice), end(x_slice));
   double min_y = *min_element(begin(y_slice), end(y_slice));
   double max_y = *max_element(begin(y_slice), end(y_slice));


   double pixel_size = .05 ;
   vector<double> xp, yp;

   for(double i=min_x; i<=max_x ; i+=pixel_size){
    xp.push_back(i);
    //cout << i << endl;
   }

   for(double i=min_y; i<=max_y; i+=pixel_size){
    yp.push_back(i);
    //cout << i << endl;
   }

   //vector<vector<int>> votes_count;
   int votes_count[(xp.size()-1)][(yp.size()-1)] = {{0}};

    //cout << xp.size() << '\n' << yp.size() << '\n' <<  sizeof(votes_count) << '\n' << sizeof(votes_count[0]);

   for(int i=0; i < x_slice.size(); i++){
    double xi = x_slice[i] , yi = y_slice[i];

    int xin=0;
    for(int j=0; j < (xp.size()-1); j++){
        if(xi >= xp[j] && xi <= xp[j+1]){
            xin = j;
            break;
        }
    }

    int yin=0;
    for(int k=0; k < (yp.size()-1); k++){
        if(yi >= yp[k] && yi <= yp[k+1]){
            yin = k;
            break;
        }
    }

    votes_count[xin][yin] += 1;
    //cout << i << ' ' << xin << ' ' << yin << ' ' << votes_count[xin][yin] << endl;

   }

   ofstream myfile;
   myfile.open ("example.txt");
   for(int i = 0; i < xp.size()-1; i++){
        string row;
        for(int j = 0; j < yp.size()-1; j++){
           myfile << votes_count[i][j] << ' ';
        }
        myfile << '\n';
   }
  myfile.close();


    return 0 ;
}
