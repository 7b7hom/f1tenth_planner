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
#include <set>
#include <queue>
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
    int node_idx;
    int layer_idx;
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

// 스플라인 결과를 담기 위한 구조체: x, y 방향 계수, 행렬 M, 정규화된 노멀 벡터
struct SplineResult {
    MatrixXd coeffs_x;              // 각 구간의 x 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd coeffs_y;              // 각 구간의 y 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd M;                     // 스플라인 계수 계산에 사용된 시스템 행렬
    MatrixXd normvec_normalized;    // 각 구간의 법선 벡터를 정규화한 값 (구간 개수 x 2)
};

typedef vector<double> DVector;
typedef vector<int>    IVector;
typedef map<string, DVector> DMap;
typedef map<string, IVector> IMap;
typedef vector<vector<Node>> NodeMap;
typedef tuple<int, int> ITuple;

typedef pair<int, int> IPair; // <layerIdx, nodeIdx>
typedef vector<IPair> IPairVector; // 엣지 연결 여부 확인용 value vector
typedef map<IPair, IPairVector> IPairAdjList; // key: 기준 노드, value: key와 연결된 다음 레이어의 노드 인덱스 IPair
typedef map<IPair, map<IPair, Spline>> SplineMap;

extern DMap gtpl_map;
extern DMap sampling_map;

struct ActionSet {
    string action_id; // "straight"
    vector<MatrixXd> coeffs; // x_coeff, y_coeff
    vector<MatrixXd> path_param; // path, psi, kappa, el_lengths 
    NodeMap nodes; // [[None, None], start_node]
    vector<IPair> node_idx; // [0, path.size()-1]
};

class Graph {
private:
    bool isDirected;
public:
    IPairAdjList adjLists;
    Graph(bool directed = true);
    void addEdge(IPair srcIdx, IPair dstIdx);
    void printGraph();
    bool getChildNodes(const IPair& parentIdx, IPairVector& childIdx);
    bool getParentNodes(const IPair& childIdx, IPairVector& parentIdx, int num_layers);
    void removeEdge(const IPair& srcIdx, const IPair& dstIdx, SplineMap* splineMap, int& remove_cnt, int num_layers);
};

// visualization.cpp
void plotHeading(const DVector &x, const DVector &y, const DVector &psi, double scale);
void plotHeading(const NodeMap& nodesPerLayer, double scale);
void visual(const NodeMap& nodesPerLayer, Graph& graph, const Offline_Params& params);

// helper_func.cpp
double normalizeAngle(double angle);
void readDMapFromCSV(const string& pathname, DMap& map);
void writeDMapToCSV(const string& pathname, DMap& map, char delimiter = ',');
void map_size(DMap& map);
void addDVectorToMap(DMap& map, string attr, const IVector* idx_array = nullptr);
void samplePointsFromRaceline(const DVector& kappa,     // 곡률
                              const DVector& dist,      // 점 사이 거리
                              double d_curve,           // 곡선 구간 샘플링 간격
                              double d_straight,        // 직선 구간 샘플링 간격
                              double curve_th,          // 곡선 판단 기준 곡률
                              IVector& idx_array);

//genSplines.cpp
void calcHeading(DVector &x_raceline,
                 DVector &y_raceline, 
                 DVector &psi);

VectorXd computeEuclideanDistances(const MatrixXd& path);
SplineResult calcSplines(const Node& startNode, const Node& endNode);

void calcCurvature(NodeMap& nodesPerLayer);
bool checkKappaValidity(const Vector4d& coeffs_x,
                        const Vector4d& coeffs_y,
                        const VectorXd& t_steps,
                        double max_allowed_kappa);
Vector2d computeSplinePosition(const RowVector4d& coeff_x, const RowVector4d& coeff_y, double t);
                        
void genNode(NodeMap& nodesPerLayer,        
            IVector& raceline_index_array,  
            const double veh_width,
            float lat_resolution);
void genEdge(Graph& graph, 
    SplineMap &splineMap,
    const NodeMap& nodesPerLayer, 
    const Offline_Params& params,
    const IVector& raceline_index_array,
    bool closed);