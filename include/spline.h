#pragma once
// #include "graph_planner.hpp"
#include "graph.hpp"

unique_ptr<Spline> calcSplines(const MatrixXd &path,
                               double psi_s,
                               double psi_e,
                               bool use_dist_scaling = true);
pair<Graph, SplineMap> genEdges(DMap &stMap,
                                NodeMap &nodesPerLayer,
                                const IVector &raceline_index_array,
                                YAML::Node &params);
                                
pair<vector<Vector2d>, VectorXd> samplingSpline(MatrixXd &coeffs_x, MatrixXd &coeffs_y, YAML::Node &params);

void pruneEdge(SplineMap &splineMap,
               Graph &wayptGraph,
               NodeMap &nodesPerLayer);