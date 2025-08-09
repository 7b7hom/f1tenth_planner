#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>  
#include <cmath>
#include <algorithm>
#include <time.h>
#include <memory>
#include <unordered_set>
#include <omp.h>

#include <Eigen/Dense>
#include "rapidcsv.h"
#include "matplotlibcpp.h"
#include "config_berlin.h"


//////////////////////////////////////////////////////////////////////
// Index names used in gtplMap
// #define INDEX_NAME  "COLUME_NAME_in_TRAJ_LTPL_CSV_FILE"
//////////////////////////////////////////////////////////////////////

// Index names originated from TRAJ_LTPL_CSV_FILE
#define POS_X       "x_ref_m"
#define POS_Y       " y_ref_m" // Do not delete the first white space in the following!!!
#define WIDTH_L     " width_left_m"
#define WIDTH_R     " width_right_m"
#define NORM_X      " x_normvec_m"
#define NORM_Y      " y_normvec_m"
#define NORM_L      " alpha_m"
#define RL_KAPPA    " kappa_racetraj_radpm"
#define RL_S        " s_racetraj_m"
#define RL_PSI      " psi_racetraj_rad"
#define RL_VX       " vx_racetraj_mps"

// Index names additionally created for this program
#define RL_dS       "delta_s"
#define RL_X        "x_raceline"
#define RL_Y        "y_raceline"
#define LB_X        "x_bound_l"
#define LB_Y        "y_bound_l"
#define RB_X        "x_bound_r"
#define RB_Y        "y_bound_r"
#define LB_PSI      "psi_bound_l"
#define RB_PSI      "psi_bound_r"

using namespace std;
using namespace rapidcsv;
using namespace Eigen;
namespace plt = matplotlibcpp;

struct Node {
    double x;
    double y;
    double psi;
    bool raceline;
};

struct Spline {
    MatrixXd coeffs_x;          
    MatrixXd coeffs_y;          
    VectorXd kappa;
    VectorXd el_lengths;   
    double cost;
    bool raceline;
};

typedef vector<double> DVector;
typedef vector<int>    IVector;
typedef map<string, DVector> DMap;
typedef map<string, IVector> IMap;
typedef vector<vector<Node>> NodeMap;

typedef pair<int, int> IPair; // <layerIdx, nodeIdx>
typedef vector<IPair> IPairVector; // 엣지 연결 여부 확인용 value vector
typedef map<IPair, IPairVector> IPairAdjList; // key: 기준 노드, value: key와 연결된 다음 레이어의 노드 인덱스 IPair
typedef map<IPair, map<IPair, Spline>> SplineMap;

struct ActionSet {
    string action_id; // "straight"
    vector<MatrixXd> coeffs; // x_coeff, y_coeff
    vector<MatrixXd> path_param; // path, psi, kappa, el_lengths 
    NodeMap nodes; // [[None, None], start_node]
    vector<IPair> node_idx; // [0, path.size()-1]
};

// visualization.cpp
void plotHeading(const DVector &x, const DVector &y, const DVector &psi, double scale);
void plotHeading(const NodeMap& nodesPerLayer, double scale);
void plotAllSplines(const IPairAdjList& edgeList, const SplineMap& splineMap, const string &color);
void plotSpline(const Spline& spline, const string& color);

// helper_func.cpp
unique_ptr<string> Load(const string& filename);
double normalizeAngle(double angle);
pair<DVector, DVector> computeBoundRight(DVector &pos_x, DVector &pos_y,
                                         DVector &norm_x, DVector &norm_y,
                                         DVector &width_r);
pair<DVector, DVector> computeBoundLeft(DVector &pos_x, DVector &pos_y,
                                        DVector &norm_x, DVector &norm_y,
                                        DVector &width_l);
pair<DVector, DVector> computeRaceline(DVector &pos_x, DVector &pos_y,
                                       DVector &norm_x, DVector &norm_y,
                                       DVector &norm_l);

DVector computeDeltaS(DVector &rl_s);
DMap readDMapFromCSV(const string& pathname);
void writeDMapToCSV(const string& pathname, DMap& map, char delimiter = ',');
void map_size(DMap& map);
bool checkInsideBounds(DMap &stMap, const Vector2d& pos, const float veh_width);
void printSplineInfo(const SplineMap& splineMap, const NodeMap& nodesPerLayer);
DVector calcHeading(DVector &x_raceline, DVector &y_raceline);

//genSplines.cpp

void visual(DMap &gtMap, DMap &stMap, const NodeMap &nodesPerLayer, const SplineMap &splineMap);
