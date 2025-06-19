#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>  
#include <cmath>
#include "rapidcsv.h"
#include "matplotlibcpp.h"

#define __x_ref "x_ref_m"
#define __y_ref " y_ref_m"
#define __width_right " width_right_m"
#define __width_left " width_left_m"
#define __x_normvec " x_normvec_m"
#define __y_normvec " y_normvec_m"
#define __alpha " alpha_m"
#define __kappa " kappa_racetraj_radpm"
#define __s_racetraj " s_racetraj_m"
#define __x_raceline "x_raceline"
#define __y_raceline "y_raceline"
#define __x_bound_r "x_bound_r"
#define __y_bound_r "y_bound_r"
#define __x_bound_l "x_bound_l"
#define __y_bound_l "y_bound_l"
#define __x_sampling "x_sampling"
#define __y_sampling "y_sampling"
#define __delta_s "delta_s"

using namespace std;
using namespace rapidcsv;
namespace plt = matplotlibcpp;

// 유지보수를 위하여 typedef로 선언 
typedef vector<double> DVector;
typedef vector<int>    IVector;
typedef map<string, DVector> DMap;
typedef map<string, IVector> IMap;

struct Node {
    int layer_idx;
    double x;
    double y;
    double psi;
    double kappa;
    bool raceline;
};

void readDmapFromCSV(const string& pathname, DMap& map);
void writeDMaptoCSV(const string& pathname, const DMap& map, char delimiter);

