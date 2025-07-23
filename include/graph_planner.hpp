#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>  
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include "config.h"
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
#define __vx " vx_racetraj_mps"

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

struct Node {
    double x;
    double y;
    double psi;
    double kappa;
    bool raceline;
};

// 스플라인 결과
struct Spline {
    MatrixXd coeffs_x;          // x 계수
    MatrixXd coeffs_y;          // y 계수  
    VectorXd el_lengths;   
};

typedef vector<double> DVector;
typedef vector<int>    IVector;
typedef map<string, DVector> DMap;
typedef map<string, IVector> IMap;

typedef vector<vector<Node>> NodeMap;
// typedef tuple<int, int> ITuple;
// typedef map<ITuple, IVector> TupleMap;

typedef pair<int, int> IPair; // <layerIdx, nodeIdx>
typedef vector<IPair> IPairVector; // 엣지 연결 여부 확인용 value vector
typedef map<IPair, IPairVector> IPairAdjList; // key: 기준 노드, value: key와 연결된 다음 레이어의 노드 인덱스 IPair

// typedef pair<IPair, IPair> EdgeKey;
// typedef map<EdgeKey, SplineResult> SplineMap;

// nodeMap[layerIdx][nodeIdx]
// splineMap[layerIdx][nodeIdx] = spline 
// typedef map<IPair, map<IPair, Spline, myCompare>, myCompare> SplineMap;
typedef map<IPair, map<IPair, Spline>> SplineMap;

#define LayerIdx(pair) (pair.first)
#define NodeIdx(pair) (pair.second)

class Graph {
private:
    bool isDirected;

public:
    IPairAdjList adjLists;
    Graph(bool directed = true);
    void addEdge(IPair srcIdx, IPair dstIdx);
    void printGraph();
    void getChildNodes(IPair& parentIdx, IPairVector& childIdx);
    void getParentNodes(IPair& childIdx, IPairVector& parentIdx);
    void removeEdge(IPair& srcIdx, IPair& dstIdx);
};

extern DMap gtpl_map;
extern DMap sampling_map;

// plotting
void plotHeading(const DVector &x, const DVector &y, const DVector &psi, double scale);
void plotHeading(const NodeMap& nodesPerLayer, double scale);
void plotAllSplines(const IPairAdjList& edgeList, const SplineMap& splineMap, const string &color);
void visual(const Graph& edgeList, const NodeMap& nodesPerLayer, const SplineMap& splineMap, const string &color);

// toCSV
void readDMapFromCSV(const string& pathname, DMap& map);
void writeDMapToCSV(const string& pathname, DMap& map, char delimiter = ',');
void map_size(DMap& map);

//genSplines
unique_ptr<Spline> calcSplines(const MatrixXd &path,
                                     double psi_s, 
                                     double psi_e, 
                                     bool use_dist_scaling);
void genEdges(NodeMap &nodesPerLayer, 
              Graph &edgeList,
              SplineMap &splineMap,
              const IVector &raceline_index_array,
              const float lat_offset,
              const float lat_resolution,
              const float curve_thr,
              const int max_lat_steps,
              const float stepsize_approx,
              const float min_vel_race,
              const float max_lateral_accel,
              const float veh_turn);