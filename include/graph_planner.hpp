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
#include <set>
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

struct Node {
    int layer_idx;
    int node_idx;
    double x;
    double y;
    double psi;
    double kappa;
    bool raceline;
};

struct EdgeInfo {
    double offline_cost;
    double spline_len;
    Eigen::RowVector4d coeffs_x_orig; // 이 엣지를 생성한 스플라인의 X 계수
    Eigen::RowVector4d coeffs_y_orig; // 이 엣지를 생성한 스플라인의 Y 계수

    // 기본 생성자 (필드 초기화)
    EdgeInfo() : offline_cost(0.0), spline_len(0.0) {}
};

struct SplineResult {
    MatrixXd coeffs_x;
    MatrixXd coeffs_y;
    MatrixXd M;
    MatrixXd normvec_normalized;
    VectorXd ds;
};

struct SplinePoint {
    double x, y;
    double psi;
    double kappa;
};

typedef vector<double> DVector; // Double type의 벡터
typedef vector<int>    IVector; // Int type의 벡터
typedef map<string, DVector> DMap;
typedef map<string, IVector> IMap;
typedef vector<vector<Node>> NodeMap;
typedef tuple<int, int> ITuple;
typedef map<ITuple, std::map<int, EdgeInfo>> GraphAdjListsMap;

class Graph {
private:
    GraphAdjListsMap adjLists;
    bool isDirected;

public:
    Graph(bool directed = true);
    void addEdge(ITuple srcKey, int destIdx);
    void addEdge(ITuple srcKey, int destIdx, const Eigen::RowVector4d& coeffs_x, const Eigen::RowVector4d& coeffs_y, double spline_len);
    void printGraph();
    void getChildIdx(ITuple srcKey, IVector& childIdx);
    void getParentNode(int target_layer, int value, vector<ITuple>& parent, int num_layers);
    void removeEdge(const ITuple& parent, int value); 
    const GraphAdjListsMap& getAdjLists() const;
    GraphAdjListsMap& getAdjLists_mutable(); 
};

