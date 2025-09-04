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
#include <utility>
#include "rapidcsv.h"
#include "matplotlibcpp.h"
#include "config_modena.h"
#include <map>
#include <tuple>

#define __x_ref "x_ref_m"
#define __y_ref "y_ref_m"
#define __width_right "width_right_m"
#define __width_left "width_left_m"
#define __x_normvec "x_normvec_m"
#define __y_normvec "y_normvec_m"
#define __alpha "alpha_m"
#define __kappa "kappa_racetraj_radpm"
#define __s_racetraj "s_racetraj_m"
#define __psi "psi_racetraj_rad"

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

typedef vector<double> DVector; // Double type의 벡터
typedef vector<int>    IVector; // Int type의 벡터
typedef map<string, DVector> DMap;
typedef map<string, IVector> IMap;

typedef vector<vector<Node>> NodeMap;
typedef tuple<int, int> ITuple;
typedef map<ITuple, std::map<int, EdgeInfo>> GraphAdjListsMap;


using XYPair = pair<DVector, DVector>;

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

struct StartPosResult {
    ITuple start_key{ITuple(-1,-1)};
    int    next_idx{-1};
    bool   out_of_track{true};
};

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

struct LayerParams{
    int layer_idx;
    double wL, wR, alpha;
    Vector2d ref_xy;
    Vector2d norm_vec;

    int raceline_index;
    double start_alpha;
    int num_nodes;
};

// -------- function prototypes --------
// I/O, 전처리
DMap    readDMapFromCSV(const string& pathname);
size_t  writeDMapToCSV(const string& pathname, DMap& map, char delimiter);
XYPair  computeBoundRight(const DMap& src);
XYPair  computeBoundLeft(const DMap& src);
XYPair  computeBoundRace(const DMap& src);
DVector computeDeltaS(const DVector& s, bool closed);
IVector samplePointsFromRaceline(const DVector& kappa, const DVector& dist,
                                 double d_curve, double d_straight, double curve_th);
double  normalizeAngle(double angle);
DVector calcHeading(const DVector& x, const DVector& y);

// 노드/스플라인/검증
vector<LayerParams> computeLayerParams(const DMap& sampled_map, double veh_width, float lat_resolution);
NodeMap buildNodeGrid(const vector<LayerParams>& layers, const DMap& sampled_map, float lat_resolution);
NodeMap fillNodeHeadings(NodeMap& nodes, const vector<LayerParams>& layers, const DMap& map);

SplineResult calcSplines(const MatrixXd& path,
                         const VectorXd* el_lengths_ptr,
                         double psi_s, double psi_e, bool use_dist_scaling);
SplinePoint  evaluateSpline(const RowVector4d& coeff_x, const RowVector4d& coeff_y,
                            double t, double ds_current, bool normalized_t);
bool isPointInsideTrackBounds(const DMap& sampled_map, double x, double y);
bool checkSplineValidity(const RowVector4d& coeff_x, const RowVector4d& coeff_y,
                         double ds_current, const Offline_Params& params,
                         const DMap& sampled_map);

// 그래프 생성/정리/코스트
Graph makeEmptyGraph(const NodeMap& nodesPerLayer);
Graph addRacelineEdges(Graph& graph, const NodeMap& nodes, const Offline_Params& params, const DMap& sampled_map);
Graph addCandidateEdges(Graph& graph, const NodeMap& nodes, const Offline_Params& params, const DMap& sampled_map);

Graph prune_graph(Graph graph, int num_layers, bool closed);
Graph gen_offline_cost(Graph graph, const Offline_Params& params, const NodeMap& nodesPerLayer);

// 시작점/필수 엣지/시각화
StartPosResult set_startpos(const Vector2d& pos_est, double heading_est,
                            const Offline_Params& params, const NodeMap& nodesPerLayer,
                            const DMap& sampled_map, double max_heading_offset_rad);
Graph ensure_edge_between(Graph graph, const NodeMap& nodesPerLayer, const Offline_Params& params,
                          const ITuple& start_key, int next_idx, const DMap& sampled_map);
void visual(const DMap& gtpl_map, const DMap& sampled_map,
            const NodeMap& nodesPerLayer, Graph& graph, const Offline_Params& params,
            const Vector2d& pos_est, double heading_est, const ITuple& start_key, int next_idx);

void dump_available_keys(const DMap& m);
void assert_has_keys(const DMap& m, const std::vector<std::string>& keys, const char* where);
