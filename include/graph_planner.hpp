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
#include "config_modena.h"


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

// extern DMap gtpl_map;
extern DMap stMap;

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
void plotAllSplines(const IPairAdjList& edgeList, const SplineMap& splineMap, const string &color);
void plotSpline(const Spline& spline, const string& color);
void visual(DMap &gtpl_map,
            const NodeMap &nodesPerLayer,
            const SplineMap &splineMap);

// helper_func.cpp
unique_ptr<string> Load(const string& filename);
double normalizeAngle(double angle);
void readDMapFromCSV(const string& pathname, DMap& map);
void writeDMapToCSV(const string& pathname, DMap& map, char delimiter = ',');
void map_size(DMap& map);
void addDVectorToMap(DMap& map, string attr);
bool checkInsideBounds(const Vector2d& pos, const float veh_width);
void printSplineInfo(const SplineMap& splineMap, const NodeMap& nodesPerLayer);

//genSplines.cpp
void calcHeading(DVector &x_raceline,
                 DVector &y_raceline, 
                 DVector &psi);
unique_ptr<Spline> calcSplines(const MatrixXd &path,
                                     double psi_s, 
                                     double psi_e, 
                                     bool use_dist_scaling=true);
VectorXd calcKappa(MatrixXd &coeffs_x,
                   MatrixXd &coeffs_y,
                   VectorXd &t_steps);
pair<Graph, SplineMap> genEdges(NodeMap &nodesPerLayer, 
              const IVector &raceline_index_array,
              Offline_Params& params);

pair<VectorXd, VectorXd> interpSplines(MatrixXd &coeffs_x,
                        MatrixXd &coeffs_y,
                        float stepsize_approx,
                        const float& veh_width,
                        double spline_len = NAN,
                        int no_interp_points = 10);

struct SplineTask {
    IPair start;
    IPair end;
    MatrixXd path;
    double psi_start;
    double psi_end;
    bool is_raceline;
};