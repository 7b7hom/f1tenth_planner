#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>  
#include <cmath>
#include <memory>
#include <Eigen/Dense>
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
#define __psi " psi_racetraj_rad"

#define __x_raceline "x_raceline"
#define __y_raceline "y_raceline"
#define __x_bound_r "x_bound_r"
#define __y_bound_r "y_bound_r"
#define __x_bound_l "x_bound_l"
#define __y_bound_l "y_bound_l"
#define __delta_s "delta_s"
#define __psi_bound_l "psi_bound_l"
#define __psi_bound_r "psi_bound_r"

using namespace std;
using namespace rapidcsv;
using namespace Eigen;
namespace plt = matplotlibcpp;

typedef vector<double> DVector;
typedef vector<int>    IVector;
typedef map<string, DVector> DMap;
typedef map<string, IVector> IMap;

struct Node {
    int layer_idx;
    int node_idx;
    double x;
    double y;
    double psi;
    double kappa;
    bool raceline;
};

struct SplineResult {
    MatrixXd coeffs_x;
    MatrixXd coeffs_y;
    MatrixXd M;
    MatrixXd normvec_normalized;
};

void readDmapFromCSV(const string& pathname, DMap& map);
void writeDMaptoCSV(const string& pathname, const DMap& map, char delimiter);