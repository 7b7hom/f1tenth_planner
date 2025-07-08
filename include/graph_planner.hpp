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
#define LayerIdx(pair) (pair.first)
#define NodeIdx(pair) (pair.second)

class Graph {
private:
    IPairAdjList adjLists;
    bool isDirected;

public:
    Graph(bool directed = true);
    void addEdge(IPair srcNodeIdx, IPair destNodeIdx);
    void printGraph();
    void getChildIdx(IPair srcNodeIdx, IPairVector& childNodeIdx);
    void getParentNode(int target_layer, int value, IPairVector& parent);
    void removeEdge(IPair& parent, int value);
};
