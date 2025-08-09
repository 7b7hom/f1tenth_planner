#pragma once
#include "graph_planner.hpp"

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

unique_ptr<Spline> calcSplines(const MatrixXd &path,
                                     double psi_s, 
                                     double psi_e, 
                                     bool use_dist_scaling=true);
VectorXd calcKappa(MatrixXd &coeffs_x,
                   MatrixXd &coeffs_y,
                   VectorXd &t_steps);
pair<Graph, SplineMap> genEdges(DMap &stMap,
                                NodeMap &nodesPerLayer,
                                const IVector &raceline_index_array,
                                Offline_Params &params);

pair<VectorXd, VectorXd> interpSplines(DMap &stMap,
                                       MatrixXd &coeffs_x,
                                       MatrixXd &coeffs_y,
                                       float stepsize_approx,
                                       const float &veh_width,
                                       double spline_len = NAN,
                                       int no_interp_points = 10);